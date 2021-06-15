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
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "hifir4m_config.hpp"
// avoid sorting!
#include "HIF.hpp"
#include "hifir4m_num_array.hpp"

namespace hifir4m {
// structure of preconditioner and solver

#ifdef HIFIR4M_USE_32INT
typedef int integer_type;
#else
typedef mwSignedIndex integer_type;
#endif

class MexNspFilter : public hif::NspFilter {
 public:
  using base = hif::NspFilter;

  explicit MexNspFilter(const std::size_t start = 0,
                        const std::size_t end = static_cast<std::size_t>(-1))
      : base(start, end) {}

  virtual ~MexNspFilter() {
    if (x_in) mxDestroyArray(x_in);
  }

  std::string f_name = "#unknown#";  ///< function to evaluate
  mxArray *f_handle = nullptr;       ///< function handle
  mutable mxArray *x_in = nullptr;   ///< input "wrapper" of x

  virtual void user_filter(void *x, const std::size_t n,
                           const char dtype) const override {
    if (dtype != 'd')
      mexErrMsgIdAndTxt("hifir4m:nspFilter:badDtype",
                        "only double precision is supported");
    if (f_name == "#unknown#")
      mexErrMsgIdAndTxt("hifir4m:nspFilter:unknownFunc", "unset function name");
    if (f_name == "feval" && !f_handle)
      mexErrMsgIdAndTxt("hifir4m:nspFilter:missingHandler",
                        "missing function handle");
    if (!x_in) {
      if (dtype == 'd')
        x_in = mxCreateDoubleMatrix(n, 1, mxREAL);
      else
        x_in = mxCreateDoubleMatrix(n, 1, mxCOMPLEX);
    }
    if (n != mxGetM(x_in))
      mexErrMsgIdAndTxt("hifir4m:nspFilter:badShape", "unmatched sizes");
    // copy x to buffer
    if (dtype == 'd')
      std::copy_n(reinterpret_cast<const double *>(x), n, mxGetPr(x_in));
    else
      std::copy_n(reinterpret_cast<const std::complex<double> *>(x), n,
                  reinterpret_cast<std::complex<double> *>(mxGetData(x_in)));
    mxArray *rhs[2];
    const mwSize nrhs = f_handle ? 2 : 1;
    if (f_handle) {
      // if using function handle, then the first one is the function handle,
      // while the input argument is given following it
      rhs[0] = f_handle;
      rhs[1] = x_in;
    } else
      rhs[0] = x_in;
    // NOTE: We have to create a new array each time to ensure the safety of
    // MATLAB runtime system.
    // TODO: can we avoid this?
    mxArray *lhs;
    // call MATLAB directly
    mexCallMATLAB(1, &lhs, nrhs, rhs, f_name.c_str());
    // copy back to data to HIFIR
    if (dtype == 'd')
      std::copy_n(mxGetPr(lhs), n, reinterpret_cast<double *>(x));
    else
      std::copy_n(
          reinterpret_cast<const std::complex<double> *>(mxGetData(lhs)), n,
          reinterpret_cast<std::complex<double> *>(x));
    mxDestroyArray(lhs);  // free the array
  }

  // function to enable user override function option
  inline void enable_or() { _type = USER_OR; }
};

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
  ///< preconditioner type
  using ksp_factory_t = hif::ksp::KSPFactory<prec_t, _value_type>;
  ///< KSP factory type
  using solver_t = typename ksp_factory_t::fgmres;  ///< FMGRES type
  ///< updated operator type
  static constexpr bool IS_MIXED = IsMixed;  ///< mixed flag
  static constexpr bool IS_REAL = _IS_REAL;  ///< real flag

  // attributes
  std::shared_ptr<prec_t> M;      ///< preconditioner attribute
  std::shared_ptr<solver_t> ksp;  ///< KSP solver
};

enum {
  HIFIR4M_CREATE = 0,   ///< create database
  HIFIR4M_GET = 1,      ///< get database
  HIFIR4M_CLEAR = 2,    ///< clear database
  HIFIR4M_CHECK = 3,    ///< check emptyness
  HIFIR4M_DESTROY = 4,  ///< destroy database
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
    mexErrMsgIdAndTxt("hifir4m:options:badInput", "input must be structure");
  hif::Options opt = hif::get_default_options();
  // WARNING! The structure from MATLAB should be initialized by a structured
  // routine with the following fields that are accessed. In addition, the
  // index must align with the order in C++ struct of Options
  auto get_field = [&](const int field) -> double {
    return mxGetScalar(mxGetFieldByNumber(rhs, 0, field));
  };
  opt.tau_L = get_field(0);
  opt.tau_U = get_field(1);
  opt.kappa_d = get_field(2);
  opt.kappa = get_field(3);
  opt.alpha_L = get_field(4);
  opt.alpha_U = get_field(5);
  opt.rho = get_field(6);
  opt.c_d = get_field(7);
  opt.c_h = get_field(8);
  opt.N = get_field(9);
  opt.verbose = get_field(10);
  opt.rf_par = get_field(11);
  opt.reorder = get_field(12);
  opt.spd = get_field(13);
  opt.check = get_field(14);
  opt.pre_scale = get_field(15);
  opt.symm_pre_lvls = get_field(16);
  opt.threads = get_field(17);
  opt.mumps_blr = get_field(18);
  opt.fat_schur_1st = get_field(19);
  opt.rrqr_cond = get_field(20);
  opt.pivot = get_field(21);
  opt.gamma = get_field(22);
  opt.beta = get_field(23);
  opt.is_symm = get_field(24);
  opt.no_pre = get_field(25);
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
  if (*rptr == 1) {
    // convert to zero based
    for (mwSize i = 0; i < n + 1; ++i) --rptr[i];
    for (mwSize i = 0; i < (mwSize)nnz; ++i) --cptr[i];
  }
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
  const auto elim_nz = A.eliminate(1e-15, true);
  if (opts.verbose) hif_info("eliminated %zd small entries in A", elim_nz);
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

#if 0
/**
 * @brief Accessing inv(M) with 2 RHS (for complex conjugate pair)
 *
 * @tparam IsMixed Wether or not the database uses mixed precision
 * @tparam ValueType Data type used, e.g., \a double or \a complex<double>
 * @param[in] id ID tag of the database
 * @param[in] rhs right-hand side vector (nx2)
 * @param[out] lhs left-hand side result, i.e., lhs=inv(A)*rhs (nx2)
 * @param[in] r (optional) rank for the final Schur complement
 * @return double overhead-free wall-clock time
 */
template <bool IsMixed, class ValueType = double>
inline double M_solve2(int id, const mxArray *rhs, mxArray *lhs,
                       const std::size_t r = 0u) {
  auto data = database<IsMixed, ValueType>(HIFIR4M_GET, id);

  if (!data->M)
    mexErrMsgIdAndTxt("hifir4m:M_solve2:emptyM", "M has not yet factorized");

  if (mxGetM(lhs) != 2u || mxGetM(rhs) != 2u)
    mexErrMsgIdAndTxt("hifir4m:M_solve2:badRhsSize", "rhs/lhs must be 2-by-n");
  if (mxGetN(rhs) != mxGetN(lhs) || mxGetN(rhs) != data->M->nrows())
    mexErrMsgIdAndTxt("hifir4m:M_solve2:badRhsSize",
                      "rhs size does not agree with lhs or M");

  using array_t = hif::Array<std::array<ValueType, 2>>;
  array_t b(data->M->nrows(), (std::array<ValueType, 2> *)mxGetData(rhs), true),
      x(data->M->nrows(), (std::array<ValueType, 2> *)mxGetData(lhs), true);
  hif::DefaultTimer timer;
  timer.start();
  try {
    data->M->solve_mrhs(b, x, r);
  } catch (const std::exception &e) {
    mexErrMsgIdAndTxt("hifir4m:M_solve:failedSolve",
                      "M_solve failed with message:\n%s", e.what());
  }
  timer.finish();
  return timer.time();  // give M solve time to the user
}
#endif

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

// KSP solve

/**
 * @brief Solving with FGMRES solver
 *
 * @tparam IsMixed Wether or not the database uses mixed precision
 * @tparam ValueType Data type used, e.g., \a double or \a complex<double>
 * @param[in] id ID tag of the database
 * @param[in] restart Restart of GMRES (30)
 * @param[in] max_iter Maximum iteration allowed (500)
 * @param[in] rtol Relative tolerance for residual convergence (1e-6)
 * @param[in] verbose Verbose flag (true)
 * @param[in] rowptr rowptr row pointer mex array (int32)
 * @param[in] colind colind column index mex array (int32)
 * @param[in] val value mex array (double)
 * @param[in] rhs Right-hand side b vector
 * @param[in,out] lhs Left-hand side solution and initial guess vector
 * @return A tuple of flag, iterations, final residual and overhead-free
 *          wall-clock time are returned in this routine.
 */
template <bool IsMixed, class ValueType = double>
inline std::tuple<int, int, double, double> KSP_solve(
    int id, const int restart, const int max_iter, const double rtol,
    const bool verbose, const mxArray *rowptr, const mxArray *colind,
    const mxArray *val, const mxArray *rhs, mxArray *lhs,
    const int iter_refines = 1, const int *cst_nsp = nullptr,
    const char *fname = nullptr, const mxArray *fhdl = nullptr) {
  using ksp_t = typename HIFIR4M_Database<IsMixed, ValueType>::solver_t;
  if (!std::is_floating_point<ValueType>::value && !mxIsComplex(rhs) &&
      !mxIsComplex(lhs) && !mxIsComplex(val))
    mexErrMsgIdAndTxt("hifir4m:M_solve:complex", "input array is not complex.");

  auto data = database<IsMixed, ValueType>(HIFIR4M_GET, id);
  if (!data->M)
    mexErrMsgIdAndTxt("hifir4m:KSP_solve:emptyM", "M has not yet factorized");

  if (mxGetN(lhs) != 1u || mxGetN(rhs) != 1u)
    mexErrMsgIdAndTxt("hifir4m:KSP_solve:badRhsSize",
                      "rhs/lhs must be column vector");
  // get matrix
  mwSize n;
  integer_type *rptr, *cptr;
  ValueType *vptr;
  std::vector<ValueType> buf;
  convert_crs_mx2pointer(std::string("hifir4m:KSP_solve"), rowptr, colind,
                         &rptr, &cptr, &n);
  get_num_data(val, &vptr, buf);

  if (mxGetM(rhs) != mxGetM(lhs) || mxGetM(rhs) != data->M->nrows() ||
      data->M->nrows() != n)
    mexErrMsgIdAndTxt("hifir4m:KSP_solve:badRhsSize",
                      "rhs size does not agree with lhs, M or A");

  data->ksp.reset(new ksp_t(data->M));
  auto &ksp = *data->ksp;
  ksp.set_M(data->M);  // setup preconditioner

  // create csr wrapper from HIFIR
  hif::CRS<ValueType, integer_type> A(n, n, rptr, cptr, vptr, true);
  // arrays
  using array_t = hif::Array<ValueType>;
  ValueType *rhs_ptr, *lhs_ptr;
  // only needed in non-interleaved complex
  std::vector<ValueType> rhs_buf, lhs_buf;
  get_num_data(rhs, &rhs_ptr, rhs_buf);
  get_num_data(lhs, &lhs_ptr, lhs_buf);
  array_t b(n, rhs_ptr, true), x(n, lhs_ptr, true);

  if (ksp.is_arnoldi() && restart > 0) ksp.set_restart_or_cycle(restart);
  if (max_iter > 0) ksp.set_maxit(max_iter);
  if (rtol > 0.0) ksp.set_rtol(rtol);
  const int irs = iter_refines > 0 ? iter_refines : 1;
  ksp.set_inner_steps(irs);

  // enable const null space filter
  if (cst_nsp)
    data->M->nsp.reset(new MexNspFilter(cst_nsp[0] - 1, cst_nsp[1]));
  else if (fname) {
    if (verbose)
      hif_info(
          "setting up user (right) null space filter with\n"
          "\tcallback name: %s\n"
          "\tfunction handle: %s",
          fname, (fhdl ? "yes" : "no"));
    // setup user filter
    data->M->nsp.reset(new MexNspFilter());
    MexNspFilter *mex_nsp = dynamic_cast<MexNspFilter *>(data->M->nsp.get());
    if (!mex_nsp)
      mexErrMsgIdAndTxt("hifir4m:KSP_solve:badNsp",
                        "failed to dynamic_cast nullspace filter");
    mex_nsp->f_name = fname;
    if (fhdl) mex_nsp->f_handle = const_cast<mxArray *>(fhdl);
    mex_nsp->enable_or();
  }

  hif::DefaultTimer timer;
  int flag;
  std::size_t iters;
  const int kernel_type =
      irs <= 1 ? hif::ksp::TRADITION : hif::ksp::ITERATIVE_REFINE;
  timer.start();
  try {
    std::tie(flag, iters) =
        ksp.solve(A, b, x, kernel_type, true /* always with guess */, verbose);
  } catch (const std::exception &e) {
    mexErrMsgIdAndTxt("hifir4m:KSP_solve:failedSolve",
                      "KSP_solve failed with message:\n%s", e.what());
  }
  timer.finish();
  const double tt = timer.time();
  const double res = data->ksp->get_resids().back();
  data->ksp.reset();                       // free
  if (data->M->nsp) data->M->nsp.reset();  // release const nullspace filter
  set_num_data(lhs, lhs_ptr);
  return std::make_tuple(flag, (int)iters, res, tt);
}

#if 0
/**
 * @brief Solve left null space with GMRES
 *
 * @tparam IsMixed Wether or not the database uses mixed precision
 * @tparam ValueType Data type used, e.g., \a double or \a complex<double>
 * @tparam UseHi using hi-precision kernels
 * @param[in] id ID tag of the database
 * @param[in] restart Restart of GMRES (30)
 * @param[in] max_iter Maximum iteration allowed (500)
 * @param[in] rtol Relative tolerance for residual convergence (1e-6)
 * @param[in] verbose Verbose flag (true)
 * @param[in] rowptr rowptr row pointer mex array (int32)
 * @param[in] colind colind column index mex array (int32)
 * @param[in] val value mex array (double)
 * @param[in] rhs Right-hand side b vector
 * @param[in,out] lhs Left-hand side solution and initial guess vector
 * @return A tuple of flag, iterations, final residual and overhead-free
 *          wall-clock time are returned in this routine.
 */
template <bool IsMixed, bool UseHi, class ValueType = double>
inline std::tuple<int, int, double, double> KSP_null_solve(
    int id, const int restart, const int max_iter, const double rtol,
    const bool verbose, const mxArray *rowptr, const mxArray *colind,
    const mxArray *val, const mxArray *rhs, mxArray *lhs) {
  using ksp_t = typename HIFIR4M_Database<IsMixed, ValueType>::null_solver_t;
  using ksp_hi_t =
      typename HIFIR4M_Database<IsMixed, ValueType>::null_hi_solver_t;

  auto data = database<IsMixed, ValueType>(HIFIR4M_GET, id);
  if (!data->M)
    mexErrMsgIdAndTxt("hifir4m:KSP_null_solve:emptyM",
                      "M has not yet factorized");

  if (mxGetN(lhs) != 1u || mxGetN(rhs) != 1u)
    mexErrMsgIdAndTxt("hifir4m:KSP_null_solve:badRhsSize",
                      "rhs/lhs must be column vector");
  // get matrix
  mwSize n;
  integer_type *rptr, *cptr;
  ValueType *vptr;
  convert_crs_mx2pointer(std::string("hifir4m:KSP_solve"), rowptr, colind, val,
                         &rptr, &cptr, &vptr, &n);

  if (mxGetM(rhs) != mxGetM(lhs) || mxGetM(rhs) != data->M->nrows() ||
      data->M->nrows() != n)
    mexErrMsgIdAndTxt("hifir4m:KSP_solve:badRhsSize",
                      "rhs size does not agree with lhs, M or A");

  // create csr wrapper from HIFIR
  hif::CRS<ValueType, integer_type> A(n, n, rptr, cptr, vptr, true);
  // arrays
  using array_t = hif::Array<ValueType>;
  array_t b(n, (ValueType *)mxGetData(rhs), true),
      x(n, (ValueType *)mxGetData(lhs), true);

  if (!UseHi) {
    data->ksp_null.reset(new ksp_t(data->M));
    auto &ksp = *data->ksp_null;
    ksp.set_M(data->M);  // setup preconditioner
    if (ksp.is_arnoldi() && restart > 0) ksp.set_restart_or_cycle(restart);
    if (max_iter > 0) ksp.set_maxit(max_iter);
    if (rtol > 0.0) ksp.set_rtol(rtol);
  } else {
    data->ksp_null_hi.reset(new ksp_hi_t(data->M));
    auto &ksp = *data->ksp_null_hi;
    ksp.set_M(data->M);  // setup preconditioner
    if (ksp.is_arnoldi() && restart > 0) ksp.set_restart_or_cycle(restart);
    if (max_iter > 0) ksp.set_maxit(max_iter);
    if (rtol > 0.0) ksp.set_rtol(rtol);
  }

  hif::DefaultTimer timer;
  int flag;
  std::size_t iters;
  timer.start();
  try {
    if (!UseHi)
      std::tie(flag, iters) = data->ksp_null->solve(
          A, b, x, hif::ksp::TRADITION, true /* always with guess */, verbose);
    else
      std::tie(flag, iters) = data->ksp_null_hi->solve(
          A, b, x, hif::ksp::TRADITION, true /* always with guess */, verbose);
  } catch (const std::exception &e) {
    mexErrMsgIdAndTxt("hifir4m:KSP_solve:failedSolve",
                      "KSP_solve failed with message:\n%s", e.what());
  }
  timer.finish();
  const double tt = timer.time();
  const double res = !UseHi ? data->ksp_null->get_resids().back()
                            : data->ksp_null_hi->get_resids().back();
  !UseHi ? data->ksp_null.reset() : data->ksp_null_hi.reset();  // free
  return std::make_tuple(flag, (int)iters, res, tt);
}

#endif

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

  if (!data->M->can_export()) {
    mexWarnMsgIdAndTxt(
        "hifir4m:M_export:cannot_export",
        "Exporting data is not allowed. 1) compiled with sparse last level, or "
        "2) using interval-based data stuctures for L and U");
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
