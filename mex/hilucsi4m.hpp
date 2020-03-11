/*
                This file is part of HILUCSI4M project

    Copyright (C) 2019 NumGeom Group at Stony Brook University

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef HILUCSI4M_HPP_
#define HILUCSI4M_HPP_

#include <algorithm>
#include <exception>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "hilucsi4m_config.hpp"
// avoid sorting!
#include "HILUCSI.hpp"

namespace hilucsi4m {
// structure of preconditioner and solver

template <class MType, class ValueType = void>
class OperatorUpdateSolver
    : public hilucsi::ksp::UserOperatorBase<MType, ValueType> {
  using _base = hilucsi::ksp::UserOperatorBase<MType, ValueType>;

 public:
  using array_type = typename _base::array_type;
  using M_type = typename _base::M_type;
  using size_type = typename _base::size_type;

  OperatorUpdateSolver() = delete;

  explicit OperatorUpdateSolver(const M_type &M) : _base(M), _v(M.nrows()) {}

  using _base::solve;

  template <class Matrix>
  inline void solve(const Matrix &A, const array_type &b, const size_type,
                    array_type &x0) const {
    // step 1: (2*I-A*inv(M)) * b
    _base::solve(b, _v);
    hilucsi::mt::mv_nt(A, _v, x0);
    const size_type n = _v.size();
    for (size_type i(0); i < n; ++i) _v[i] = 2.0 * b[i] - x0[i];
    // step 2: inv(M)*v
    _base::solve(_v, x0);
  }

 protected:
  mutable array_type _v;  ///< workspace
};

class MexNspFilter : public hilucsi::NspFilter {
 public:
  using base = hilucsi::NspFilter;

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
      mexErrMsgIdAndTxt("hilucsi4m:nspFilter:badDtype",
                        "only double precision is supported");
    if (f_name == "#unknown#")
      mexErrMsgIdAndTxt("hilucsi4m:nspFilter:unknownFunc",
                        "unset function name");
    if (f_name == "feval" && !f_handle)
      mexErrMsgIdAndTxt("hilucsi4m:nspFilter:missingHandler",
                        "missing function handle");
    if (!x_in) x_in = mxCreateDoubleMatrix(n, 1, mxREAL);
    if (n != mxGetM(x_in))
      mexErrMsgIdAndTxt("hilucsi4m:nspFilter:badShape", "unmatched sizes");
    // copy x to buffer
    std::copy_n(reinterpret_cast<const double *>(x), n, mxGetPr(x_in));
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
    // copy back to data to HILUCSI
    std::copy_n(mxGetPr(lhs), n, reinterpret_cast<double *>(x));
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
template <bool IsMixed>
struct HILUCSI4M_Database {
 private:
  using _value_type = typename std::conditional<IsMixed, double, void>::type;

 public:
  using prec_t =
      typename std::conditional<IsMixed, hilucsi::HILUCSI<float, int>,
                                hilucsi::DefaultHILUCSI>::type;
  ///< preconditioner type
  using ksp_factory_t = hilucsi::ksp::KSPFactory<prec_t, _value_type>;
  ///< KSP factory type
  using solver_t = typename ksp_factory_t::fgmres;           ///< FMGRES type
  using null_solver_t = typename ksp_factory_t::gmres_null;  ///< left null
  using update_operator_t = OperatorUpdateSolver<prec_t, _value_type>;
  ///< updated operator type
  static constexpr bool IS_MIXED = IsMixed;  ///< mixed flag

  // attributes
  std::shared_ptr<prec_t> M;                ///< preconditioner attribute
  std::shared_ptr<solver_t> ksp;            ///< KSP solver
  std::shared_ptr<null_solver_t> ksp_null;  ///< KSP null solver
};

// // structure of preconditioner with mixed precision
// template <>
// struct HILUCSI4M_Database<true> {
//   using prec_t        = hilucsi::HILUCSI<float, int>;  ///< preconditioner
//   type using ksp_factory_t = hilucsi::ksp::KSPFactory<prec_t, double>;
//   ///< KSP factory type
//   using solver_t = ksp_factory_t::fgmres;
//   ///< FMGRES type
//   using update_operator_t = OperatorUpdateSolver<prec_t, double>;

//   std::shared_ptr<prec_t>   M;    ///< preconditioner attribute
//   std::shared_ptr<solver_t> ksp;  ///< KSP solver
// };

enum {
  HILUCSI4M_CREATE = 0,   ///< create database
  HILUCSI4M_GET = 1,      ///< get database
  HILUCSI4M_CLEAR = 2,    ///< clear database
  HILUCSI4M_DESTROY = 3,  ///< destroy database
};

/**
 * @brief unified interface for database query/creation/deletion
 *
 * @tparam IsMixed Wether or not the database uses mixed precision
 * @param[in] action Action flag
 * @param[in,out] id If action is creation, then output the new ID
 * @return A pointer of to the database is returned
 */
template <bool IsMixed>
inline HILUCSI4M_Database<IsMixed> *database(const int action, int &id) {
  // create database
  using data_t = HILUCSI4M_Database<IsMixed>;
  static std::vector<std::shared_ptr<data_t>> pool;

  switch (action) {
    case HILUCSI4M_CREATE: {
      id = pool.size();
      pool.emplace_back(new data_t);
      return pool.back().get();
    } break;
    case HILUCSI4M_GET: {
      // query database
      if (id < 0 || id >= (int)pool.size())
        mexErrMsgIdAndTxt("hilucsi4m:database:getInvalidID",
                          "%d is an invalid database ID", id);
      if (!pool[id])
        mexErrMsgIdAndTxt("hilucsi4m:database:getDeadID",
                          "%d database has been destroyed", id);
      return pool[id].get();
    } break;
    case HILUCSI4M_CLEAR: {
      if (id < 0 || id >= (int)pool.size())
        mexErrMsgIdAndTxt("hilucsi4m:database:getInvalidID",
                          "%d is an invalid database ID", id);
      if (!pool[id])
        mexErrMsgIdAndTxt("hilucsi4m:database:getDeadID",
                          "%d database has been destroyed", id);
      pool[id]->M->clear();
      return pool[id].get();
    } break;
    default: {
      if (action != HILUCSI4M_DESTROY)
        mexErrMsgIdAndTxt("hilucsi4m:database:badAction",
                          "%d is an invalid action for database", action);
      if (id < 0 || id >= (int)pool.size())
        mexErrMsgIdAndTxt("hilucsi4m:database:destroyInvalidID",
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
 * @return hilucsi::Options control parameters for HILUCSI
 *
 * @note This routine requires the field ordering in the matlab structure
 *       follows the sames order as how the attributes appear in the C++
 *       struct definition.
 */
inline hilucsi::Options create_opt_from_struct(const mxArray *rhs) {
  if (!mxIsStruct(rhs))
    mexErrMsgIdAndTxt("hilucsi4m:options:badInput", "input must be structure");
  hilucsi::Options opt = hilucsi::get_default_options();
  // WARNING! The structure from MATLAB should be initialized by a structured
  // routine with the following fields that are accessed. In addition, the
  // index must align with the order in C++ struct of Options
  auto get_field = [&](const int field) -> double {
    return mxGetScalar(mxGetFieldByNumber(rhs, 0, field));
  };
  opt.tau_L = get_field(0);
  opt.tau_U = get_field(1);
  opt.tau_d = get_field(2);
  opt.tau_kappa = get_field(3);
  opt.alpha_L = get_field(4);
  opt.alpha_U = get_field(5);
  opt.rho = get_field(6);
  opt.c_d = get_field(7);
  opt.c_h = get_field(8);
  opt.N = get_field(9);
  opt.verbose = get_field(10);
  opt.rf_par = get_field(11);
  opt.reorder = get_field(12);
  opt.saddle = get_field(13);
  opt.check = get_field(14);
  opt.pre_scale = get_field(15);
  opt.symm_pre_lvls = get_field(16);
  opt.threads = get_field(17);
  opt.fat_schur_1st = get_field(18);
  return opt;
}

// helper function to parse CSR matrix

/**
 * @brief Helper function to convert MATLAB arrays of CRS to raw pointer
 *
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
                                   const mxArray *val, int **rptr_, int **cptr_,
                                   double **vptr_, mwSize *n_) {
  if (!mxIsInt32(rowptr) || !mxIsInt32(colind))
    mexErrMsgIdAndTxt((prefix + ":badIntDtype").c_str(),
                      "rowptr and colind must be int32");
  const auto n = mxGetM(rowptr) - 1;
  int *rptr = (int *)mxGetData(rowptr), *cptr = (int *)mxGetData(colind);
  const auto nnz = rptr[n] - rptr[0];
  if (mxGetM(colind) < (mwSize)nnz || mxGetM(val) < (mwSize)nnz)
    mexErrMsgIdAndTxt((prefix + ":badLength").c_str(), "bad nnz length %d",
                      nnz);
  if (*rptr == 1) {
    // convert to zero based
    for (mwSize i = 0; i < n + 1; ++i) --rptr[i];
    for (int i = 0; i < nnz; ++i) --cptr[i];
  }
  *rptr_ = rptr;
  *cptr_ = cptr;
  *vptr_ = (double *)mxGetPr(val);
  *n_ = n;
}

// factorization

/**
 * @brief Factorize an HILUCSI instance
 *
 * @tparam IsMixed Wether or not the database uses mixed precision
 * @param[in] id ID tag of the database
 * @param[in] rowptr row pointer mex array (int32)
 * @param[in] colind column index mex array (int32)
 * @param[in] val value mex array (double)
 * @param[in] opt control parameter (1x1 structure)
 * @return overhead-free factorization wall-clock time
 */
template <bool IsMixed>
inline double factorize(int id, const mxArray *rowptr, const mxArray *colind,
                        const mxArray *val, const mxArray *opt) {
  using data_t = HILUCSI4M_Database<IsMixed>;

  auto data = database<IsMixed>(HILUCSI4M_GET, id);

  mwSize n;
  int *rptr, *cptr;
  double *vptr;
  convert_crs_mx2pointer(std::string("hilucsi4m:factorize"), rowptr, colind,
                         val, &rptr, &cptr, &vptr, &n);

  // create csr wrapper from HILUCSI
  hilucsi::CRS<double, int> A(n, n, rptr, cptr, vptr, true);
  // get options
  const auto opts = create_opt_from_struct(opt);
  if (!data->M) data->M.reset(new typename data_t::prec_t());
  hilucsi::DefaultTimer timer;
  timer.start();
  try {
    data->M->factorize(A, 0, opts);
  } catch (const std::exception &e) {
    mexErrMsgIdAndTxt("hilucsi4m:factorize:failedFac",
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
 * @param[in] id ID tag of the database
 * @param[in] rhs right-hand side vector
 * @param[out] lhs left-hand side result, i.e., lhs=inv(A)*rhs
 * @return double overhead-free wall-clock time
 */
template <bool IsMixed>
inline double M_solve(int id, const mxArray *rhs, mxArray *lhs) {
  auto data = database<IsMixed>(HILUCSI4M_GET, id);

  if (!data->M)
    mexErrMsgIdAndTxt("hilucsi4m:M_solve:emptyM", "M has not yet factorized");

  if (mxGetN(lhs) != 1u || mxGetN(rhs) != 1u)
    mexErrMsgIdAndTxt("hilucsi4m:M_solve:badRhsSize",
                      "rhs/lhs must be column vector");
  if (mxGetM(rhs) != mxGetM(lhs) || mxGetM(rhs) != data->M->nrows())
    mexErrMsgIdAndTxt("hilucsi4m:M_solve:badRhsSize",
                      "rhs size does not agree with lhs or M");

  using array_t = hilucsi::Array<double>;
  array_t b(data->M->nrows(), mxGetPr(rhs), true),
      x(data->M->nrows(), mxGetPr(lhs));
  hilucsi::DefaultTimer timer;
  timer.start();
  try {
    data->M->solve(b, x);
  } catch (const std::exception &e) {
    mexErrMsgIdAndTxt("hilucsi4m:M_solve:failedSolve",
                      "M_solve failed with message:\n%s", e.what());
  }
  timer.finish();
  return timer.time();  // give M solve time to the user
}

// KSP solve

/**
 * @brief Solving with FGMRES solver
 *
 * @tparam IsMixed Wether or not the database uses mixed precision
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
template <bool IsMixed>
inline std::tuple<int, int, double, double> KSP_solve(
    int id, const int restart, const int max_iter, const double rtol,
    const bool verbose, const mxArray *rowptr, const mxArray *colind,
    const mxArray *val, const mxArray *rhs, mxArray *lhs,
    const bool update = false, const int *cst_nsp = nullptr,
    const char *fname = nullptr, const mxArray *fhdl = nullptr) {
  using ksp_t = typename HILUCSI4M_Database<IsMixed>::solver_t;

  auto data = database<IsMixed>(HILUCSI4M_GET, id);
  if (!data->M)
    mexErrMsgIdAndTxt("hilucsi4m:KSP_solve:emptyM", "M has not yet factorized");

  if (mxGetN(lhs) != 1u || mxGetN(rhs) != 1u)
    mexErrMsgIdAndTxt("hilucsi4m:KSP_solve:badRhsSize",
                      "rhs/lhs must be column vector");
  // get matrix
  mwSize n;
  int *rptr, *cptr;
  double *vptr;
  convert_crs_mx2pointer(std::string("hilucsi4m:KSP_solve"), rowptr, colind,
                         val, &rptr, &cptr, &vptr, &n);

  if (mxGetM(rhs) != mxGetM(lhs) || mxGetM(rhs) != data->M->nrows() ||
      data->M->nrows() != n)
    mexErrMsgIdAndTxt("hilucsi4m:KSP_solve:badRhsSize",
                      "rhs size does not agree with lhs, M or A");

  data->ksp.reset(new ksp_t(data->M));
  auto &ksp = *data->ksp;
  ksp.set_M(data->M);  // setup preconditioner

  // create csr wrapper from HILUCSI
  hilucsi::CRS<double, int> A(n, n, rptr, cptr, vptr, true);
  // arrays
  using array_t = hilucsi::Array<double>;
  array_t b(n, mxGetPr(rhs), true), x(n, mxGetPr(lhs), true);

  if (ksp.is_arnoldi() && restart > 0) ksp.set_restart_or_cycle(restart);
  if (max_iter > 0) ksp.set_maxit(max_iter);
  if (rtol > 0.0) ksp.set_rtol(rtol);

  // enable const null space filter
  if (cst_nsp)
    data->M->nsp.reset(new MexNspFilter(cst_nsp[0] - 1, cst_nsp[1]));
  else if (fname) {
    if (verbose)
      hilucsi_info(
          "setting up user (right) null space filter with\n"
          "\tcallback name: %s\n"
          "\tfunction handle: %s",
          fname, (fhdl ? "yes" : "no"));
    // setup user filter
    data->M->nsp.reset(new MexNspFilter());
    MexNspFilter *mex_nsp = dynamic_cast<MexNspFilter *>(data->M->nsp.get());
    if (!mex_nsp)
      mexErrMsgIdAndTxt("hilucsi4m:KSP_solve:badNsp",
                        "failed to dynamic_cast nullspace filter");
    mex_nsp->f_name = fname;
    if (fhdl) mex_nsp->f_handle = const_cast<mxArray *>(fhdl);
    mex_nsp->enable_or();
  }

  hilucsi::DefaultTimer timer;
  int flag;
  std::size_t iters;
  timer.start();
  try {
    if (!update)
      std::tie(flag, iters) = ksp.solve(A, b, x, hilucsi::ksp::TRADITION,
                                        true /* always with guess */, verbose);
    else {
      const auto M =
          typename HILUCSI4M_Database<IsMixed>::update_operator_t(*data->M);
      std::tie(flag, iters) = ksp.solve_user(A, M, b, x, true, verbose);
    }
  } catch (const std::exception &e) {
    mexErrMsgIdAndTxt("hilucsi4m:KSP_solve:failedSolve",
                      "KSP_solve failed with message:\n%s", e.what());
  }
  timer.finish();
  const double tt = timer.time();
  const double res = data->ksp->get_resids().back();
  data->ksp.reset();                       // free
  if (data->M->nsp) data->M->nsp.reset();  // release const nullspace filter
  return std::make_tuple(flag, (int)iters, res, tt);
}

/**
 * @brief Solve left null space with GMRES
 *
 * @tparam IsMixed Wether or not the database uses mixed precision
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
template <bool IsMixed>
inline std::tuple<int, int, double, double> KSP_null_solve(
    int id, const int restart, const int max_iter, const double rtol,
    const bool verbose, const mxArray *rowptr, const mxArray *colind,
    const mxArray *val, const mxArray *rhs, mxArray *lhs) {
  using ksp_t = typename HILUCSI4M_Database<IsMixed>::null_solver_t;

  auto data = database<IsMixed>(HILUCSI4M_GET, id);
  if (!data->M)
    mexErrMsgIdAndTxt("hilucsi4m:KSP_null_solve:emptyM",
                      "M has not yet factorized");

  if (mxGetN(lhs) != 1u || mxGetN(rhs) != 1u)
    mexErrMsgIdAndTxt("hilucsi4m:KSP_null_solve:badRhsSize",
                      "rhs/lhs must be column vector");
  // get matrix
  mwSize n;
  int *rptr, *cptr;
  double *vptr;
  convert_crs_mx2pointer(std::string("hilucsi4m:KSP_solve"), rowptr, colind,
                         val, &rptr, &cptr, &vptr, &n);

  if (mxGetM(rhs) != mxGetM(lhs) || mxGetM(rhs) != data->M->nrows() ||
      data->M->nrows() != n)
    mexErrMsgIdAndTxt("hilucsi4m:KSP_solve:badRhsSize",
                      "rhs size does not agree with lhs, M or A");

  data->ksp_null.reset(new ksp_t(data->M));
  auto &ksp = *data->ksp_null;
  ksp.set_M(data->M);  // setup preconditioner

  // create csr wrapper from HILUCSI
  hilucsi::CRS<double, int> A(n, n, rptr, cptr, vptr, true);
  // arrays
  using array_t = hilucsi::Array<double>;
  array_t b(n, mxGetPr(rhs), true), x(n, mxGetPr(lhs), true);

  if (ksp.is_arnoldi() && restart > 0) ksp.set_restart_or_cycle(restart);
  if (max_iter > 0) ksp.set_maxit(max_iter);
  if (rtol > 0.0) ksp.set_rtol(rtol);

  hilucsi::DefaultTimer timer;
  int flag;
  std::size_t iters;
  timer.start();
  try {
    std::tie(flag, iters) = ksp.solve(A, b, x, hilucsi::ksp::TRADITION,
                                      true /* always with guess */, verbose);
  } catch (const std::exception &e) {
    mexErrMsgIdAndTxt("hilucsi4m:KSP_solve:failedSolve",
                      "KSP_solve failed with message:\n%s", e.what());
  }
  timer.finish();
  const double tt = timer.time();
  const double res = data->ksp_null->get_resids().back();
  data->ksp_null.reset();  // free
  return std::make_tuple(flag, (int)iters, res, tt);
}
}  // namespace hilucsi4m

#endif  // HILUCSI4M_HPP_
