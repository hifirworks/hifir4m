/*
                This file is part of HIFIR4M project

    Copyright (C) 2019--2021 NumGeom Group at Stony Brook University

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef HIFIR4M_HPP_
#define HIFIR4M_HPP_

#include <algorithm>
#include <array>
#include <complex>
#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "hifir4m_config.hpp"
#include "hifir4m_num_array.hpp"

namespace hifir4m {
// structure of preconditioner and solver

#ifdef HIFIR4M_USE_32INT
typedef int integer_type;
#else
typedef INT64_T integer_type;
#endif

/**
 * @brief Database structure
 *
 * @tparam IsMixed Flag indicate if or not the database uses mixed precision
 */
template <bool IsMixed, class ValueType = double>
struct HIFIR4M_Database {
 private:
  using _value_type = typename std::conditional<IsMixed, ValueType, void>::type;
  static constexpr bool _IS_REAL = std::is_floating_point<ValueType>::value;
  ///< real flag
  using _reduce_type = typename std::conditional<
      IsMixed,
      typename std::conditional<_IS_REAL, float, std::complex<float>>::type,
      typename std::conditional<_IS_REAL, double,
                                std::complex<double>>::type>::type;

 public:
  using prec_t =
      typename std::conditional<IsMixed, hif::HIF<_reduce_type, integer_type>,
                                hif::HIF<ValueType, integer_type>>::type;
  ///< updated operator type
  static constexpr bool IS_MIXED = IsMixed;  ///< mixed flag
  static constexpr bool IS_REAL = _IS_REAL;  ///< real flag

  // attributes
  std::shared_ptr<prec_t> M;  ///< preconditioner attribute
};

enum {
  HIFIR4M_CREATE = 0,       ///< create database
  HIFIR4M_GET = 1,          ///< get database
  HIFIR4M_CLEAR = 2,        ///< clear database
  HIFIR4M_CHECK = 3,        ///< check emptyness
  HIFIR4M_DESTROY = 4,      ///< destroy database
  HIFIR4M_FACTORIZE = 5,    ///< Factorize
  HIFIR4M_M_SOLVE = 6,      ///< preconditioner solve
  HIFIR4M_M_MULTIPLY = 7,   ///< matrix vector product
  HIFIR4M_EXPORT_DATA = 8,  ///< export data
  HIFIR4M_QUERY = 9,        ///< query factorization info
};

/**
 * @brief unified interface for database query/creation/deletion
 *
 * @tparam IsMixed Wether or not the database uses mixed precision
 * @tparam ValueType Data value type, e.g., \a double or \a complex<double>
 * @param[in] action Action flag
 * @param[in,out] id If action is creation, then output the new ID
 * @return A pointer of to the database is returned
 */
template <bool IsMixed, class ValueType = double>
inline HIFIR4M_Database<IsMixed, ValueType> *database(const int action,
                                                      int &id) {
  // create database
  using data_t = HIFIR4M_Database<IsMixed, ValueType>;
  static std::vector<std::shared_ptr<data_t>> pool;

  switch (action) {
    case HIFIR4M_CREATE: {
      id = pool.size();
      pool.emplace_back(new data_t);
      return pool.back().get();
    } break;
    case HIFIR4M_GET: {
      // query database
      if (id < 0 || id >= (int)pool.size())
        mexErrMsgIdAndTxt("hifir4m:database:getInvalidID",
                          "%d is an invalid database ID", id);
      if (!pool[id]) pool[id] = std::make_shared<data_t>();
      return pool[id].get();
    } break;
    case HIFIR4M_CLEAR: {
      if (id < 0 || id >= (int)pool.size())
        mexErrMsgIdAndTxt("hifir4m:database:getInvalidID",
                          "%d is an invalid database ID", id);
      if (!pool[id])
        mexErrMsgIdAndTxt("hifir4m:database:getDeadID",
                          "%d database has been destroyed", id);
      if (pool[id]->M) pool[id]->M->clear();
      return pool[id].get();
    } break;
    default: {
      if (action != HIFIR4M_DESTROY)
        mexErrMsgIdAndTxt("hifir4m:database:badAction",
                          "%d is an invalid action for database", action);
      if (id < 0 || id >= (int)pool.size())
        mexErrMsgIdAndTxt("hifir4m:database:destroyInvalidID",
                          "%d is an invalid database ID", id);
      pool[id].reset();
      return nullptr;
    } break;
  }
}

/**
 * @brief Create a opt from a MATLAB structure
 *
 * @param[in] rhs MATLAB mex array
 * @return hif::Options control parameters for HIFIR
 *
 * @note This routine requires the field ordering in the matlab structure
 *       follows the sames order as how the attributes appear in the C++
 *       struct definition.
 */
inline hif::Options create_opt_from_struct(const mxArray *rhs) {
  if (!mxIsStruct(rhs))
    mexErrMsgIdAndTxt("hifir4m:params:badInput", "input must be structure");
  hif::Options opt = hif::get_default_options();
  // WARNING! The structure from MATLAB should be initialized by a structured
  // routine with the following fields that are accessed. In addition, the
  // index must align with the order in C++ struct of Options
  auto get_field = [&](const char *field) -> double {
    auto *fd = mxGetField(rhs, 0, field);
    if (!fd)
      mexErrMsgIdAndTxt("hifir4m:params:invalidField", "%s is unknown field",
                        field);
    return mxGetScalar(fd);
  };
  opt.tau_L = get_field("tau_L");
  opt.tau_U = get_field("tau_U");
  opt.kappa_d = get_field("kappa_d");
  opt.kappa = get_field("kappa");
  opt.alpha_L = get_field("alpha_L");
  opt.alpha_U = get_field("alpha_U");
  opt.rho = get_field("rho");
  opt.c_d = get_field("c_d");
  opt.c_h = get_field("c_h");
  opt.N = get_field("N");
  opt.verbose = get_field("verbose");
  opt.rf_par = get_field("rf_par");
  opt.reorder = get_field("reorder");
  opt.spd = get_field("spd");
  opt.check = get_field("check");
  opt.pre_scale = get_field("pre_scale");
  opt.symm_pre_lvls = get_field("symm_pre_lvls");
  opt.threads = get_field("threads");
  opt.mumps_blr = get_field("mumps_blr");
  opt.fat_schur_1st = get_field("fat_schur_1st");
  opt.rrqr_cond = get_field("rrqr_cond");
  opt.pivot = get_field("pivot");
  opt.gamma = get_field("gamma");
  opt.beta = get_field("beta");
  opt.is_symm = get_field("is_symm");
  opt.no_pre = get_field("no_pre");
  return opt;
}

// helper function to parse CSR matrix

/**
 * @brief Helper function to convert MATLAB arrays of CRS to raw pointer
 *
 * @tparam ValueType Data type used, e.g., \a double or \a complex<double>
 * @param[in] prefix MATLAB message ID prefix (for error output purpose)
 * @param[in] rowptr row pointer mex array
 * @param[in] colind column index mex array
 * @param[in] val value mex array
 * @param[out] rptr_ output pointer to the \a rowptr array
 * @param[out] cptr_ output pointer to the \a column array
 * @param[out] vptr_ output pointer to the \a val array
 * @param[out] n_ size of the system
 *
 * @note This routine assumes 0-based index of rowptr/column, in addition, the
 *       mex arrays must maintain valid while accessing the output pointers.
 */
inline void convert_crs_mx2pointer(const std::string &prefix,
                                   const mxArray *rowptr, const mxArray *colind,
                                   integer_type **rptr_, integer_type **cptr_,
                                   mwSize *n_) {
#ifdef HILIUCSI4M_USE_32INT
  if (!mxIsInt32(rowptr) || !mxIsInt32(colind))
    mexErrMsgIdAndTxt((prefix + ":badIntDtype").c_str(),
                      "rowptr and colind must be int32");
#else
  if (!mxIsInt64(rowptr) || !mxIsInt64(colind))
    mexErrMsgIdAndTxt((prefix + ":badIntDtype").c_str(),
                      "rowptr and colind must be int64");
#endif
  const auto n = mxGetM(rowptr) - 1;
  integer_type *rptr = (integer_type *)mxGetData(rowptr),
               *cptr = (integer_type *)mxGetData(colind);
  const auto nnz = rptr[n] - rptr[0];
  if (mxGetM(colind) < (mwSize)nnz)
    mexErrMsgIdAndTxt((prefix + ":badLength").c_str(), "bad nnz length %d",
                      nnz);
  if (*rptr != 1 && *rptr != 0)
    mexErrMsgIdAndTxt((prefix + ":badIndexBase").c_str(),
                      "must be {0,1}-based");
  *rptr_ = rptr;
  *cptr_ = cptr;
  *n_ = n;
}

// factorization

/**
 * @brief Factorize an HIFIR instance
 *
 * @tparam IsMixed Wether or not the database uses mixed precision
 * @tparam ValueType Data type used, e.g., \a double or \a complex<double>
 * @param[in] id ID tag of the database
 * @param[in] rowptr row pointer mex array (int32)
 * @param[in] colind column index mex array (int32)
 * @param[in] val value mex array (double)
 * @param[in] opt control parameter (1x1 structure)
 * @return overhead-free factorization wall-clock time
 */
template <bool IsMixed, class ValueType = double>
inline double factorize(int id, const mxArray *rowptr, const mxArray *colind,
                        const mxArray *val, const mxArray *opt) {
  using data_t = HIFIR4M_Database<IsMixed, ValueType>;

  if (!std::is_floating_point<ValueType>::value && !mxIsComplex(val))
    mexErrMsgIdAndTxt("hifir4m:factorize:complex",
                      "input val array is not complex.");

  auto data = database<IsMixed, ValueType>(HIFIR4M_GET, id);

  mwSize n;
  integer_type *rptr, *cptr;
  ValueType *vptr;
  std::vector<ValueType> buf;  // only used in complex and no interleaved
  convert_crs_mx2pointer(std::string("hifir4m:factorize"), rowptr, colind,
                         &rptr, &cptr, &n);
  get_num_data(val, &vptr, buf);

  // create csr wrapper from HIFIR
  hif::CRS<ValueType, integer_type> A(n, n, rptr, cptr, vptr, true);
  // get options
  const auto opts = create_opt_from_struct(opt);
  if (!data->M) data->M.reset(new typename data_t::prec_t());
  hif::DefaultTimer timer;
  timer.start();
  try {
    data->M->factorize(A, opts);
  } catch (const std::exception &e) {
    data->M->clear();
    mexErrMsgIdAndTxt("hifir4m:factorize:failedFac",
                      "factorization failed with message:\n%s", e.what());
  }
  timer.finish();
  return timer.time();  // give factorization time to the user
}

// M solver (preconditioner solve)

/**
 * @brief Accessing inv(M)
 *
 * @tparam IsMixed Wether or not the database uses mixed precision
 * @tparam ValueType Data type used, e.g., \a double or \a complex<double>
 * @param[in] id ID tag of the database
 * @param[in] rhs right-hand side vector
 * @param[out] lhs left-hand side result, i.e., lhs=inv(A)*rhs
 * @param[in] trans (optional) transpose/Hermitian flag, default is \a false
 * @param[in] r (optional) rank for the final Schur complement
 * @return double overhead-free wall-clock time
 */
template <bool IsMixed, class ValueType = double>
inline double M_solve(int id, const mxArray *rhs, mxArray *lhs,
                      const bool trans = false, const std::size_t r = 0u) {
  auto data = database<IsMixed, ValueType>(HIFIR4M_GET, id);

  if (!data->M)
    mexErrMsgIdAndTxt("hifir4m:M_solve:emptyM", "M has not yet factorized");
  if (!std::is_floating_point<ValueType>::value && !mxIsComplex(rhs) &&
      !mxIsComplex(lhs))
    mexErrMsgIdAndTxt("hifir4m:M_solve:complex", "input array is not complex.");

  if (mxGetN(lhs) != 1u || mxGetN(rhs) != 1u)
    mexErrMsgIdAndTxt("hifir4m:M_solve:badRhsSize",
                      "rhs/lhs must be column vector");
  if (mxGetM(rhs) != mxGetM(lhs) || mxGetM(rhs) != data->M->nrows())
    mexErrMsgIdAndTxt("hifir4m:M_solve:badRhsSize",
                      "rhs size does not agree with lhs or M");

  using array_t = hif::Array<ValueType>;
  ValueType *rhs_ptr, *lhs_ptr;
  // only needed in non-interleaved complex
  std::vector<ValueType> rhs_buf, lhs_buf;
  get_num_data(rhs, &rhs_ptr, rhs_buf);
  get_num_data(lhs, &lhs_ptr, lhs_buf);
  array_t b(data->M->nrows(), rhs_ptr, true),
      x(data->M->nrows(), lhs_ptr, true);
  hif::DefaultTimer timer;
  timer.start();
  try {
    data->M->solve(b, x, trans, r);
  } catch (const std::exception &e) {
    mexErrMsgIdAndTxt("hifir4m:M_solve:failedSolve",
                      "M_solve failed with message:\n%s", e.what());
  }
  timer.finish();
  set_num_data(lhs, lhs_ptr);
  return timer.time();  // give M solve time to the user
}

// M solve with inner iteration

/**
 * @brief Factorize an HIFIR instance
 *
 * @tparam IsMixed Wether or not the database uses mixed precision
 * @tparam ValueType Data type used, e.g., \a double or \a complex<double>
 * @param[in] id ID tag of the database
 * @param[in] rowptr row pointer mex array (int32)
 * @param[in] colind column index mex array (int32)
 * @param[in] val value mex array (double)
 * @param[in] rhs right-hand side vector
 * @param[in] N inner iteration steps
 * @param[out] lhs left-hand side result
 * @param[in] trans (optional) transpose/Hermitian flag, default is \a false
 * @param[in] r (optional) rank for the final Schur complement
 * @return overhead-free factorization wall-clock time
 */
template <bool IsMixed, class ValueType = double>
inline double M_solve(int id, const mxArray *rowptr, const mxArray *colind,
                      const mxArray *val, const mxArray *rhs, const int N,
                      mxArray *lhs, const bool trans = false,
                      const std::size_t r = static_cast<std::size_t>(-1)) {
  if (N < 2) return M_solve<IsMixed, ValueType>(id, rhs, lhs);  // quick return

  if (!std::is_floating_point<ValueType>::value && !mxIsComplex(rhs) &&
      !mxIsComplex(lhs) && !mxIsComplex(val))
    mexErrMsgIdAndTxt("hifir4m:M_solve:complex", "input array is not complex.");

  auto data = database<IsMixed, ValueType>(HIFIR4M_GET, id);

  if (!data->M)
    mexErrMsgIdAndTxt("hifir4m:M_solve_inner:emptyM",
                      "M has not yet factorized");

  if (mxGetN(lhs) != 1u || mxGetN(rhs) != 1u)
    mexErrMsgIdAndTxt("hifir4m:M_solve_inner:badRhsSize",
                      "rhs/lhs must be column vector");
  if (mxGetM(rhs) != mxGetM(lhs) || mxGetM(rhs) != data->M->nrows())
    mexErrMsgIdAndTxt("hifir4m:M_solve_inner:badRhsSize",
                      "rhs size does not agree with lhs or M");

  mwSize n;
  integer_type *rptr, *cptr;
  ValueType *vptr;
  std::vector<ValueType> buf;
  convert_crs_mx2pointer(std::string("hifir4m:M_solve_inner"), rowptr, colind,
                         &rptr, &cptr, &n);
  get_num_data(val, &vptr, buf);

  bool local_buf = false;
  if (*rptr != 0) {
    local_buf = true;
    integer_type *ptr_bak = rptr;
    rptr = (integer_type *)mxMalloc((n + 1) * sizeof(integer_type));
    const auto index_base = ptr_bak[0];
    for (mwSize i(0); i < n + 1; ++i) rptr[i] = ptr_bak[i] - index_base;
    ptr_bak = cptr;
    const mwSize nnz = rptr[n] - index_base;
    cptr = (integer_type *)mxMalloc(sizeof(integer_type) * nnz);
    for (mwSize i(0); i < nnz; ++i) cptr[i] = ptr_bak[i] - index_base;
  }

  // create csr wrapper from HIFIR
  hif::CRS<ValueType, integer_type> A(n, n, rptr, cptr, vptr, true);
  using array_t = hif::Array<ValueType>;
  ValueType *rhs_ptr, *lhs_ptr;
  // only needed in non-interleaved complex
  std::vector<ValueType> rhs_buf, lhs_buf;
  get_num_data(rhs, &rhs_ptr, rhs_buf);
  get_num_data(lhs, &lhs_ptr, lhs_buf);
  array_t b(data->M->nrows(), rhs_ptr, true),
      x(data->M->nrows(), lhs_ptr, true);
  hif::DefaultTimer timer;
  try {
    data->M->hifir(A, b, N, x, trans, r);  // call inner iteration kernel
  } catch (const std::exception &e) {
    mexErrMsgIdAndTxt("hifir4m:M_solve_inner:failedSolve",
                      "M_solve_inner failed with message:\n%s", e.what());
  }
  timer.finish();
  set_num_data(lhs, lhs_ptr);
  if (local_buf) {
    mxFree(rptr);
    mxFree(cptr);
  }
  return timer.time();  // give M solve time to the user
}

// Multilevel matrix-vector multiplication

/**
 * @brief Perform M*x
 *
 * @tparam IsMixed Wether or not the database uses mixed precision
 * @tparam ValueType Data type used, e.g., \a double or \a complex<double>
 * @param[in] id ID tag of the database
 * @param[in] rhs right-hand side vector
 * @param[out] lhs left-hand side result, i.e., lhs=M*rhs
 * @param[in] trans (optional) transpose/Hermitian flag, default is \a false
 * @param[in] r (optional) rank for the final Schur complement
 * @return double overhead-free wall-clock time
 */
template <bool IsMixed, class ValueType = double>
inline double M_multiply(int id, const mxArray *rhs, mxArray *lhs,
                         const bool trans = false, const std::size_t r = 0u) {
  auto data = database<IsMixed, ValueType>(HIFIR4M_GET, id);

  if (!data->M)
    mexErrMsgIdAndTxt("hifir4m:M_multiply:emptyM", "M has not yet factorized");
  if (!std::is_floating_point<ValueType>::value && !mxIsComplex(rhs) &&
      !mxIsComplex(lhs))
    mexErrMsgIdAndTxt("hifir4m:M_solve:complex", "input array is not complex.");
  if (mxGetN(lhs) != 1u || mxGetN(rhs) != 1u)
    mexErrMsgIdAndTxt("hifir4m:M_multiply:badRhsSize",
                      "rhs/lhs must be column vector");
  if (mxGetM(rhs) != mxGetM(lhs) || mxGetM(rhs) != data->M->nrows())
    mexErrMsgIdAndTxt("hifir4m:M_multiply:badRhsSize",
                      "rhs size does not agree with lhs or M");

  using array_t = hif::Array<ValueType>;
  ValueType *rhs_ptr, *lhs_ptr;
  // only needed in non-interleaved complex
  std::vector<ValueType> rhs_buf, lhs_buf;
  get_num_data(rhs, &rhs_ptr, rhs_buf);
  get_num_data(lhs, &lhs_ptr, lhs_buf);
  array_t b(data->M->nrows(), rhs_ptr, true),
      x(data->M->nrows(), lhs_ptr, true);
  hif::DefaultTimer timer;
  timer.start();
  try {
    data->M->mmultiply(b, x, trans, r);
  } catch (const std::exception &e) {
    mexErrMsgIdAndTxt("hifir4m:M_multiply:failedSolve",
                      "M_multiply failed with message:\n%s", e.what());
  }
  timer.finish();
  set_num_data(lhs, lhs_ptr);
  return timer.time();  // give M solve time to the user
}

/**
 * @brief Extract the internal data from HIFIR
 *
 * @tparam IsMixed Wether or not the database uses mixed precision
 * @tparam ValueType Data type used, e.g., \a double or \a complex<double>
 * @param[in] id ID tag of the database
 * @param[in] destroy Destroy internal data while exporting them?
 * @param[out] hilu Multi-level ILU factorization
 */
template <bool IsMixed, class ValueType = double>
inline void M_export(int id, const bool destroy, mxArray **hilu) {
  using data_t = HIFIR4M_Database<IsMixed, ValueType>;
  using crs_t = typename data_t::prec_t::crs_type;
  using value_t = typename crs_t::value_type;
  static constexpr mxComplexity DTYPE =
      std::is_floating_point<ValueType>::value ? mxREAL : mxCOMPLEX;
  using scalar_type = typename hif::ValueTypeTrait<value_t>::value_type;

  auto data = database<IsMixed, ValueType>(HIFIR4M_GET, id);
  if (!data || !data->M || data->M->empty()) {
    // empty multi-level M
    *hilu = mxCreateCellMatrix(0, 0);
    return;
  }

  // attributes in crs
  static const char *crs_attrs[] = {"row_ptr", "col_ind", "val", "nrows",
                                    "ncols"};

  // attributes in a single level
  static const char *prec_attrs[] = {"L_B",   "D_B",   "U_B",   "E",
                                     "F",     "s_row", "s_col", "p",
                                     "p_inv", "q",     "q_inv"};

  value_t *vptr1, *vptr2, *vptr3, *vptr4, *vptr5;
  std::vector<value_t> vbuf1, vbuf2, vbuf3, vbuf4, vbuf5;

  // initialize cell array (hilu)
  *hilu = mxCreateCellMatrix(data->M->levels(), 1);
  int index(0);
  for (auto iter = data->M->precs().begin(); iter != data->M->precs().end();
       ++iter) {
    typename crs_t::size_type m, n, nnz_L, nnz_U, nnz_E, nnz_F;
    iter->inquire_sizes(m, n, nnz_L, nnz_U, nnz_E, nnz_F);
    // allocate for L_B
    auto *L_row_ptr =
        mxCreateUninitNumericMatrix(m + 1, 1, mxINT64_CLASS, mxREAL);
    auto *L_col_ind =
        mxCreateUninitNumericMatrix(nnz_L, 1, mxINT64_CLASS, mxREAL);
    auto *L_val = mxCreateUninitNumericMatrix(
        nnz_L, 1, IsMixed ? mxSINGLE_CLASS : mxDOUBLE_CLASS, DTYPE);
    // allocate for D_B
    auto *D_B = mxCreateUninitNumericMatrix(
        m, 1, IsMixed ? mxSINGLE_CLASS : mxDOUBLE_CLASS, DTYPE);
    // allocate for U_B
    auto *U_row_ptr =
        mxCreateUninitNumericMatrix(m + 1, 1, mxINT64_CLASS, mxREAL);
    auto *U_col_ind =
        mxCreateUninitNumericMatrix(nnz_U, 1, mxINT64_CLASS, mxREAL);
    auto *U_val = mxCreateUninitNumericMatrix(
        nnz_U, 1, IsMixed ? mxSINGLE_CLASS : mxDOUBLE_CLASS, DTYPE);
    // allocate for E
    auto *E_row_ptr =
        mxCreateUninitNumericMatrix(n - m + 1, 1, mxINT64_CLASS, mxREAL);
    auto *E_col_ind =
        mxCreateUninitNumericMatrix(nnz_E, 1, mxINT64_CLASS, mxREAL);
    auto *E_val = mxCreateUninitNumericMatrix(
        nnz_E, 1, IsMixed ? mxSINGLE_CLASS : mxDOUBLE_CLASS, DTYPE);
    // allocate for F
    auto *F_row_ptr =
        mxCreateUninitNumericMatrix(m + 1, 1, mxINT64_CLASS, mxREAL);
    auto *F_col_ind =
        mxCreateUninitNumericMatrix(nnz_F, 1, mxINT64_CLASS, mxREAL);
    auto *F_val = mxCreateUninitNumericMatrix(
        nnz_F, 1, IsMixed ? mxSINGLE_CLASS : mxDOUBLE_CLASS, DTYPE);
    // row/column scaling
    auto *s_row = mxCreateUninitNumericMatrix(
             n, 1, IsMixed ? mxSINGLE_CLASS : mxDOUBLE_CLASS, mxREAL),
         *s_col = mxCreateUninitNumericMatrix(
             n, 1, IsMixed ? mxSINGLE_CLASS : mxDOUBLE_CLASS, mxREAL);
    // row/column permutation arrays
    auto *p = mxCreateUninitNumericMatrix(n, 1, mxINT64_CLASS, mxREAL),
         *p_inv = mxCreateUninitNumericMatrix(n, 1, mxINT64_CLASS, mxREAL),
         *q = mxCreateUninitNumericMatrix(n, 1, mxINT64_CLASS, mxREAL),
         *q_inv = mxCreateUninitNumericMatrix(n, 1, mxINT64_CLASS, mxREAL);

    get_num_data(L_val, &vptr1, vbuf1);
    get_num_data(D_B, &vptr2, vbuf2);
    get_num_data(U_val, &vptr3, vbuf3);
    get_num_data(E_val, &vptr4, vbuf4);
    get_num_data(F_val, &vptr5, vbuf5);
    iter->template export_sparse_data<crs_t>(
        (integer_type *)mxGetData(L_row_ptr),
        (integer_type *)mxGetData(L_col_ind), vptr1, vptr2,
        (integer_type *)mxGetData(U_row_ptr),
        (integer_type *)mxGetData(U_col_ind), vptr3,
        (integer_type *)mxGetData(E_row_ptr),
        (integer_type *)mxGetData(E_col_ind), vptr4,
        (integer_type *)mxGetData(F_row_ptr),
        (integer_type *)mxGetData(F_col_ind), vptr5,
        (scalar_type *)mxGetData(s_row), (scalar_type *)mxGetData(s_col),
        (integer_type *)mxGetData(p), (integer_type *)mxGetData(p_inv),
        (integer_type *)mxGetData(q), (integer_type *)mxGetData(q_inv));
    set_num_data(L_val, vptr1);
    set_num_data(D_B, vptr2);
    set_num_data(U_val, vptr3);
    set_num_data(E_val, vptr4);
    set_num_data(F_val, vptr5);

    // increment index by 1 for MATLAB/Fortran index base
    std::for_each((integer_type *)mxGetData(L_row_ptr),
                  (integer_type *)mxGetData(L_row_ptr) + mxGetM(L_row_ptr),
                  [](integer_type &v) { ++v; });
    std::for_each((integer_type *)mxGetData(L_col_ind),
                  (integer_type *)mxGetData(L_col_ind) + mxGetM(L_col_ind),
                  [](integer_type &v) { ++v; });
    std::for_each((integer_type *)mxGetData(U_row_ptr),
                  (integer_type *)mxGetData(U_row_ptr) + mxGetM(U_row_ptr),
                  [](integer_type &v) { ++v; });
    std::for_each((integer_type *)mxGetData(U_col_ind),
                  (integer_type *)mxGetData(U_col_ind) + mxGetM(U_col_ind),
                  [](integer_type &v) { ++v; });
    std::for_each((integer_type *)mxGetData(E_row_ptr),
                  (integer_type *)mxGetData(E_row_ptr) + mxGetM(E_row_ptr),
                  [](integer_type &v) { ++v; });
    std::for_each((integer_type *)mxGetData(E_col_ind),
                  (integer_type *)mxGetData(E_col_ind) + mxGetM(E_col_ind),
                  [](integer_type &v) { ++v; });
    std::for_each((integer_type *)mxGetData(F_row_ptr),
                  (integer_type *)mxGetData(F_row_ptr) + mxGetM(F_row_ptr),
                  [](integer_type &v) { ++v; });
    std::for_each((integer_type *)mxGetData(F_col_ind),
                  (integer_type *)mxGetData(F_col_ind) + mxGetM(F_col_ind),
                  [](integer_type &v) { ++v; });
    std::for_each((integer_type *)mxGetData(p),
                  (integer_type *)mxGetData(p) + mxGetM(p),
                  [](integer_type &v) { ++v; });
    std::for_each((integer_type *)mxGetData(p_inv),
                  (integer_type *)mxGetData(p_inv) + mxGetM(p_inv),
                  [](integer_type &v) { ++v; });
    std::for_each((integer_type *)mxGetData(q),
                  (integer_type *)mxGetData(q) + mxGetM(q),
                  [](integer_type &v) { ++v; });
    std::for_each((integer_type *)mxGetData(q_inv),
                  (integer_type *)mxGetData(q_inv) + mxGetM(q_inv),
                  [](integer_type &v) { ++v; });

    // build CRS struct for L
    auto *L_struct = mxCreateStructMatrix(1, 1, 5, crs_attrs);
    mxSetFieldByNumber(L_struct, 0, 0, L_row_ptr);
    mxSetFieldByNumber(L_struct, 0, 1, L_col_ind);
    mxSetFieldByNumber(L_struct, 0, 2, L_val);
    do {
      auto *mm = mxCreateUninitNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
      *(integer_type *)mxGetData(mm) = m;
      mxSetFieldByNumber(L_struct, 0, 3, mm);
      auto *nn = mxCreateUninitNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
      *(integer_type *)mxGetData(nn) = m;
      mxSetFieldByNumber(L_struct, 0, 4, nn);
    } while (false);

    // build CRS struct for U
    auto *U_struct = mxCreateStructMatrix(1, 1, 5, crs_attrs);
    mxSetFieldByNumber(U_struct, 0, 0, U_row_ptr);
    mxSetFieldByNumber(U_struct, 0, 1, U_col_ind);
    mxSetFieldByNumber(U_struct, 0, 2, U_val);
    do {
      auto *mm = mxCreateUninitNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
      *(integer_type *)mxGetData(mm) = m;
      mxSetFieldByNumber(U_struct, 0, 3, mm);
      auto *nn = mxCreateUninitNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
      *(integer_type *)mxGetData(nn) = m;
      mxSetFieldByNumber(U_struct, 0, 4, nn);
    } while (false);

    // build CRS struct for E
    auto *E_struct = mxCreateStructMatrix(1, 1, 5, crs_attrs);
    mxSetFieldByNumber(E_struct, 0, 0, E_row_ptr);
    mxSetFieldByNumber(E_struct, 0, 1, E_col_ind);
    mxSetFieldByNumber(E_struct, 0, 2, E_val);
    do {
      auto *mm = mxCreateUninitNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
      *(integer_type *)mxGetData(mm) = n - m;
      mxSetFieldByNumber(E_struct, 0, 3, mm);
      auto *nn = mxCreateUninitNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
      *(integer_type *)mxGetData(nn) = m;
      mxSetFieldByNumber(E_struct, 0, 4, nn);
    } while (false);

    // build CRS struct for F
    auto *F_struct = mxCreateStructMatrix(1, 1, 5, crs_attrs);
    mxSetFieldByNumber(F_struct, 0, 0, F_row_ptr);
    mxSetFieldByNumber(F_struct, 0, 1, F_col_ind);
    mxSetFieldByNumber(F_struct, 0, 2, F_val);
    do {
      auto *mm = mxCreateUninitNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
      *(integer_type *)mxGetData(mm) = m;
      mxSetFieldByNumber(F_struct, 0, 3, mm);
      auto *nn = mxCreateUninitNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
      *(integer_type *)mxGetData(nn) = n - m;
      mxSetFieldByNumber(F_struct, 0, 4, nn);
    } while (false);

    // build prec struct
    auto *prec_struct = mxCreateStructMatrix(1, 1, 11, prec_attrs);
    mxSetFieldByNumber(prec_struct, 0, 0, L_struct);
    mxSetFieldByNumber(prec_struct, 0, 1, D_B);
    mxSetFieldByNumber(prec_struct, 0, 2, U_struct);
    mxSetFieldByNumber(prec_struct, 0, 3, E_struct);
    mxSetFieldByNumber(prec_struct, 0, 4, F_struct);
    mxSetFieldByNumber(prec_struct, 0, 5, s_row);
    mxSetFieldByNumber(prec_struct, 0, 6, s_col);
    mxSetFieldByNumber(prec_struct, 0, 7, p);
    mxSetFieldByNumber(prec_struct, 0, 8, p_inv);
    mxSetFieldByNumber(prec_struct, 0, 9, q);
    mxSetFieldByNumber(prec_struct, 0, 10, q_inv);

    // assemble to cell array for this level
    mxSetCell(*hilu, index, prec_struct);

    ++index;  // increment index
    // check for last level
    if (iter->is_last_level()) {
      typename crs_t::size_type nrows, ncols;
      value_t *dummy(nullptr);
      iter->inquire_or_export_dense(dummy, nrows, ncols, destroy);
      auto *mat = mxCreateUninitNumericMatrix(
          nrows, ncols, IsMixed ? mxSINGLE_CLASS : mxDOUBLE_CLASS, DTYPE);
      get_num_data(mat, &vptr1, vbuf1);
      iter->inquire_or_export_dense(vptr1, nrows, ncols, destroy);
      set_num_data(mat, vptr1);
      mxSetCell(*hilu, index, mat);
    }
  }
}
}  // namespace hifir4m

#endif  // HIFIR4M_HPP_
