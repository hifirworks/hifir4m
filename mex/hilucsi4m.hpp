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
  using M_type     = typename _base::M_type;
  using size_type  = typename _base::size_type;

  OperatorUpdateSolver() = delete;

  explicit OperatorUpdateSolver(const M_type &M) : _base(M), _v(M.nrows()) {}

  template <class Matrix>
  inline bool solve(const Matrix &A, const array_type &b, const size_type,
                    array_type &x0) const {
    // step 1: (2*I-A*inv(M)) * b
    _base::solve(A, b, size_type(0), _v);
    hilucsi::mt::mv_nt(A, _v, x0);
    const size_type n = _v.size();
    for (size_type i(0); i < n; ++i) _v[i] = 2.0 * b[i] - x0[i];
    // step 2: inv(M)*v
    _base::solve(A, _v, size_type(0), x0);
    return false;
  }

 protected:
  mutable array_type _v;  ///< workspace
};

/**
 * @brief Database structure
 *
 * @tparam IsMixed Flag indicate if or not the database uses mixed precision
 */
template <bool IsMixed>
struct HILUCSI4M_Database {
  using prec_t        = hilucsi::DefaultHILUCSI;  ///< preconditioner type
  using ksp_factory_t = hilucsi::ksp::KSPFactory<prec_t>;
  ///< KSP factory type
  using solver_t          = ksp_factory_t::fgmres;  ///< FMGRES type
  using update_operator_t = OperatorUpdateSolver<prec_t>;

  std::shared_ptr<prec_t> M;    ///< preconditioner attribute
  solver_t                ksp;  ///< KSP solver
};

// structure of preconditioner with mixed precision
template <>
struct HILUCSI4M_Database<true> {
  using prec_t        = hilucsi::HILUCSI<float, int>;  ///< preconditioner type
  using ksp_factory_t = hilucsi::ksp::KSPFactory<prec_t, double>;
  ///< KSP factory type
  using solver_t = ksp_factory_t::fgmres;
  ///< FMGRES type
  using update_operator_t = OperatorUpdateSolver<prec_t, double>;

  std::shared_ptr<prec_t> M;    ///< preconditioner attribute
  solver_t                ksp;  ///< KSP solver
};

enum {
  HILUCSI4M_CREATE  = 0,  ///< create database
  HILUCSI4M_GET     = 1,  ///< get database
  HILUCSI4M_DESTROY = 2,  ///< destroy database
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
  opt.tau_L         = get_field(0);
  opt.tau_U         = get_field(1);
  opt.tau_d         = get_field(2);
  opt.tau_kappa     = get_field(3);
  opt.alpha_L       = get_field(4);
  opt.alpha_U       = get_field(5);
  opt.rho           = get_field(6);
  opt.c_d           = get_field(7);
  opt.c_h           = get_field(8);
  opt.N             = get_field(9);
  opt.verbose       = get_field(10);
  opt.rf_par        = get_field(11);
  opt.reorder       = get_field(12);
  opt.saddle        = get_field(13);
  opt.check         = get_field(14);
  opt.pre_scale     = get_field(15);
  opt.symm_pre_lvls = get_field(16);
  opt.threads       = get_field(17);
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
  const auto n    = mxGetM(rowptr) - 1;
  int *      rptr = (int *)mxGetData(rowptr), *cptr = (int *)mxGetData(colind);
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
  *n_    = n;
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

  mwSize  n;
  int *   rptr, *cptr;
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
 * @return A tuple of flag, iterations and overhead-free wall-clock time are
 *         returned in this routine.
 */
template <bool IsMixed>
inline std::tuple<int, int, double> KSP_solve(
    int id, const int restart, const int max_iter, const double rtol,
    const bool verbose, const mxArray *rowptr, const mxArray *colind,
    const mxArray *val, const mxArray *rhs, mxArray *lhs,
    const bool update = false) {
  auto data = database<IsMixed>(HILUCSI4M_GET, id);
  if (!data->M)
    mexErrMsgIdAndTxt("hilucsi4m:KSP_solve:emptyM", "M has not yet factorized");

  if (mxGetN(lhs) != 1u || mxGetN(rhs) != 1u)
    mexErrMsgIdAndTxt("hilucsi4m:KSP_solve:badRhsSize",
                      "rhs/lhs must be column vector");
  // get matrix
  mwSize  n;
  int *   rptr, *cptr;
  double *vptr;
  convert_crs_mx2pointer(std::string("hilucsi4m:KSP_solve"), rowptr, colind,
                         val, &rptr, &cptr, &vptr, &n);

  if (mxGetM(rhs) != mxGetM(lhs) || mxGetM(rhs) != data->M->nrows() ||
      data->M->nrows() != n)
    mexErrMsgIdAndTxt("hilucsi4m:KSP_solve:badRhsSize",
                      "rhs size does not agree with lhs, M or A");

  auto &ksp = data->ksp;
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
  int                   flag;
  std::size_t           iters;
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
  return std::make_tuple(flag, (int)iters, tt);
}
}  // namespace hilucsi4m

#endif  // HILUCSI4M_HPP_
