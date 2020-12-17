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

#include "hilucsi4m.hpp"

namespace hilucsi4m {
enum {
  HILUCSI4M_FACTORIZE = HILUCSI4M_DESTROY + 1,     ///< Factorize
  HILUCSI4M_M_SOLVE = HILUCSI4M_FACTORIZE + 1,     ///< preconditioner solve
  HILUCSI4M_KSP_SOLVE = HILUCSI4M_M_SOLVE + 1,     ///< KSP solve
  HILUCSI4M_KSP_NULL = HILUCSI4M_KSP_SOLVE + 1,    ///< KSP for left null
  HILUCSI4M_EXPORT_DATA = HILUCSI4M_KSP_NULL + 1,  ///< export data
};
}  // namespace hilucsi4m

// structure of database fields
static const char* fnames[3] = {"id", "is_mixed", "is_real"};

// helper function for getting id and is_mixed
static bool decode_database_struct(const mxArray* rhs, int& id, bool& is_mixed,
                                   bool& is_real);

template <bool IsMixed, bool isReal>
static void get_fac_info(int id, mwSize& nnz, mwSize& nnz_ef, mwSize& deferrals,
                         mwSize& dy_deferrals, mwSize& drops,
                         mwSize& space_drops, mwSize& lvls);

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  using complex_t = std::complex<double>;

  if (nrhs < 1)
    mexErrMsgIdAndTxt("hilucsi4m:mexgateway",
                      "At least one argument must be given");
  const int action = (int)mxGetScalar(prhs[0]);
  // create
  if (action == hilucsi4m::HILUCSI4M_CREATE) {
    if (nlhs < 1)
      mexErrMsgIdAndTxt("hilucsi4m:mexgateway:create",
                        "must be at least one output");
    // determine if use mixed precision
    const int is_mixed = nrhs < 2 ? 0 : (int)mxGetScalar(prhs[1]);
    const int is_complex = nrhs < 3 ? 0 : (int)mxGetScalar(prhs[2]);
    int id;
    if (!is_mixed) {
      if (!is_complex)
        hilucsi4m::database<false, double>(hilucsi4m::HILUCSI4M_CREATE, id);
      else
        hilucsi4m::database<false, complex_t>(hilucsi4m::HILUCSI4M_CREATE, id);
    } else {
      if (!is_complex)
        hilucsi4m::database<true, double>(hilucsi4m::HILUCSI4M_CREATE, id);
      else
        hilucsi4m::database<true, complex_t>(hilucsi4m::HILUCSI4M_CREATE, id);
    }
    // create structure
    plhs[0] = mxCreateStructMatrix(1, 1, 3, fnames);
    auto id_mx = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    *(int*)mxGetData(id_mx) = id;
    mxSetFieldByNumber(plhs[0], 0, 0, id_mx);
    auto is_mixed_mx = mxCreateNumericMatrix(1, 1, mxLOGICAL_CLASS, mxREAL);
    *(mxLogical*)mxGetData(is_mixed_mx) = is_mixed;
    mxSetFieldByNumber(plhs[0], 0, 1, is_mixed_mx);
    auto is_real_mx = mxCreateNumericMatrix(1, 1, mxLOGICAL_CLASS, mxREAL);
    *(mxLogical*)mxGetData(is_real_mx) = is_complex == 0;
    mxSetFieldByNumber(plhs[0], 0, 2, is_real_mx);
    return;
  }  // end creation

  // for the rest cases, the second one must be given and be structure of
  // database created in the section above
  if (nrhs < 2)
    mexErrMsgIdAndTxt("hilucsi4m:mexgateway",
                      "the database token structure must be passed in");
  if (!mxIsStruct(prhs[1]))
    mexErrMsgIdAndTxt("hilucsi4m:mexgateway",
                      "the database token must be structure");
  int id;
  bool is_mixed, is_real;
  if (decode_database_struct(prhs[1], id, is_mixed, is_real))
    mexErrMsgIdAndTxt("hilucsi4m:mexgateway",
                      "failed to decode the database structure");
  // clear database
  if (action == hilucsi4m::HILUCSI4M_CLEAR) {
    if (!is_mixed) {
      if (is_real)
        hilucsi4m::database<false, double>(hilucsi4m::HILUCSI4M_CLEAR, id);
      else
        hilucsi4m::database<false, complex_t>(hilucsi4m::HILUCSI4M_CLEAR, id);
    } else {
      if (is_real)
        hilucsi4m::database<true, double>(hilucsi4m::HILUCSI4M_CLEAR, id);
      else
        hilucsi4m::database<true, complex_t>(hilucsi4m::HILUCSI4M_CLEAR, id);
    }
    return;
  }  // end clearing

  if (action == hilucsi4m::HILUCSI4M_CHECK) {
    if (!is_mixed) {
      if (is_real) {
        auto* dbase = hilucsi4m::database<false>(hilucsi4m::HILUCSI4M_GET, id);
        if (!dbase || !dbase->M)
          plhs[0] = mxCreateLogicalScalar(1);
        else
          plhs[0] = mxCreateLogicalScalar((mxLogical)dbase->M->empty());
      } else {
        auto* dbase =
            hilucsi4m::database<false, complex_t>(hilucsi4m::HILUCSI4M_GET, id);
        if (!dbase || !dbase->M)
          plhs[0] = mxCreateLogicalScalar(1);
        else
          plhs[0] = mxCreateLogicalScalar((mxLogical)dbase->M->empty());
      }
    } else {
      if (is_real) {
        auto* dbase = hilucsi4m::database<true>(hilucsi4m::HILUCSI4M_GET, id);
        if (!dbase || !dbase->M)
          plhs[0] = mxCreateLogicalScalar(1);
        else
          plhs[0] = mxCreateLogicalScalar((mxLogical)dbase->M->empty());
      } else {
        auto* dbase =
            hilucsi4m::database<true, complex_t>(hilucsi4m::HILUCSI4M_GET, id);
        if (!dbase || !dbase->M)
          plhs[0] = mxCreateLogicalScalar(1);
        else
          plhs[0] = mxCreateLogicalScalar((mxLogical)dbase->M->empty());
      }
    }
    return;
  }

  // destroy
  if (action == hilucsi4m::HILUCSI4M_DESTROY) {
    if (!is_mixed) {
      if (is_real)
        hilucsi4m::database<false>(hilucsi4m::HILUCSI4M_DESTROY, id);
      else
        hilucsi4m::database<false, complex_t>(hilucsi4m::HILUCSI4M_DESTROY, id);
    } else {
      if (is_real)
        hilucsi4m::database<true>(hilucsi4m::HILUCSI4M_DESTROY, id);
      else
        hilucsi4m::database<true, complex_t>(hilucsi4m::HILUCSI4M_DESTROY, id);
    }
    return;
  }  // end destroying

  // factorization
  if (action == hilucsi4m::HILUCSI4M_FACTORIZE) {
    // act, dbase, rowptr, colind, val, opt
    if (nrhs < 6)
      mexErrMsgIdAndTxt("hilucsi4m:mexgateway:fac",
                        "factorization requires 6 inputs");
    double tt;
    if (is_real)
      tt = is_mixed ? hilucsi4m::factorize<true>(id, prhs[2], prhs[3], prhs[4],
                                                 prhs[5])
                    : hilucsi4m::factorize<false>(id, prhs[2], prhs[3], prhs[4],
                                                  prhs[5]);
    else
      tt = is_mixed ? hilucsi4m::factorize<true, complex_t>(
                          id, prhs[2], prhs[3], prhs[4], prhs[5])
                    : hilucsi4m::factorize<false, complex_t>(
                          id, prhs[2], prhs[3], prhs[4], prhs[5]);
    if (nlhs) plhs[0] = mxCreateDoubleScalar(tt);
    if (nlhs > 1) {
      static const char* fac_info_names[7] = {
          "levels",       "nnz",   "nnz_EF",     "deferrals",
          "dy_deferrals", "drops", "space_drops"};
      // create an information structure
      mwSize info[7];
      if (!is_mixed) {
        if (is_real)
          get_fac_info<false, true>(id, info[1], info[2], info[3], info[4],
                                    info[5], info[6], info[0]);
        else
          get_fac_info<false, false>(id, info[1], info[2], info[3], info[4],
                                     info[5], info[6], info[0]);
      } else {
        if (is_real)
          get_fac_info<true, true>(id, info[1], info[2], info[3], info[4],
                                   info[5], info[6], info[0]);
        else
          get_fac_info<true, false>(id, info[1], info[2], info[3], info[4],
                                    info[5], info[6], info[0]);
      }
      plhs[1] = mxCreateStructMatrix(1, 1, 7, fac_info_names);
      for (int i = 0; i < 7; ++i) {
        auto mx_ptr = mxCreateUninitNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
        *(mxUint64*)mxGetData(mx_ptr) = info[i];
        mxSetFieldByNumber(plhs[1], 0, i, mx_ptr);
      }
    }
    return;
  }  // end factorization

  // M solve
  if (action == hilucsi4m::HILUCSI4M_M_SOLVE) {
    // act, dbase, b
    if (nrhs < 3)
      mexErrMsgIdAndTxt("hilucsi4m:mexgateway:m_solve",
                        "Accessing inv(M) requires 3 inputs");
    // do not use R2017b api for shared copy
    if (nlhs < 1)
      mexErrMsgIdAndTxt("hilucsi4m:mexgateway:m_solve",
                        "Accessing inv(M) requires 1 output");
    double tt;
    if (nrhs == 3) {
      plhs[0] = mxDuplicateArray(prhs[2]);
      if (is_real)
        tt = is_mixed ? hilucsi4m::M_solve<true>(id, prhs[2], plhs[0])
                      : hilucsi4m::M_solve<false>(id, prhs[2], plhs[0]);
      else
        tt = is_mixed
                 ? hilucsi4m::M_solve<true, complex_t>(id, prhs[2], plhs[0])
                 : hilucsi4m::M_solve<false, complex_t>(id, prhs[2], plhs[0]);
    } else {
      // inner iterations
      // act, dbase, b, rowptr, colind, val, N
      if (nrhs < 7)
        mexErrMsgIdAndTxt("hilucsi4m:mexgateway:m_solve_inner",
                          "Accessing inv(M) requires at least 7 inputs");
      const int N = (int)mxGetScalar(prhs[6]);
      plhs[0] = mxDuplicateArray(prhs[2]);
      if (is_real)
        tt = is_mixed ? hilucsi4m::M_solve<true>(id, prhs[3], prhs[4], prhs[5],
                                                 prhs[2], N, plhs[0])
                      : hilucsi4m::M_solve<false>(id, prhs[3], prhs[4], prhs[5],
                                                  prhs[2], N, plhs[0]);
      else
        tt = is_mixed ? hilucsi4m::M_solve<true, complex_t>(
                            id, prhs[3], prhs[4], prhs[5], prhs[2], N, plhs[0])
                      : hilucsi4m::M_solve<false, complex_t>(
                            id, prhs[3], prhs[4], prhs[5], prhs[2], N, plhs[0]);
    }
    if (nlhs > 1) plhs[1] = mxCreateDoubleScalar(tt);
    return;
  }  // end M solve

  if (action == hilucsi4m::HILUCSI4M_KSP_SOLVE) {
    // KSP solver
    // act, dbase, rowptr, colind, val, b, (restart, rtol, maxit, x0, verbose)
    if (nrhs < 6)
      mexErrMsgIdAndTxt("hilucsi4m:mexgateway:ksp_solve",
                        "KSP solver requires at least 6 inputs");
    if (nlhs < 1)
      mexErrMsgIdAndTxt("hilucsi4m:mexgateway:ksp_solve",
                        "KSP solver requires at least 1 output(s)");
    const int restart = nrhs < 7 ? 30 : (int)mxGetScalar(prhs[6]);
    const double rtol = nrhs < 8 ? 1e-6 : mxGetScalar(prhs[7]);
    const int maxit = nrhs < 9 ? 500 : (int)mxGetScalar(prhs[8]);
    if (nrhs > 9) {
      // has initial guess
      plhs[0] = mxDuplicateArray(prhs[9]);
    } else {
      plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[5]), 1,
                                     is_real ? mxREAL : mxCOMPLEX);
    }
    const bool verbose = nrhs < 11 ? true : (bool)mxGetScalar(prhs[10]);
    const bool update = nrhs < 12 ? false : (bool)mxGetScalar(prhs[11]);
    int flag, iters;
    double res, tt;
    if (nrhs < 13) {
      if (is_real) {
        if (is_mixed)
          std::tie(flag, iters, res, tt) = hilucsi4m::KSP_solve<true>(
              id, restart, maxit, rtol, verbose, prhs[2], prhs[3], prhs[4],
              prhs[5], plhs[0], update);
        else
          std::tie(flag, iters, res, tt) = hilucsi4m::KSP_solve<false>(
              id, restart, maxit, rtol, verbose, prhs[2], prhs[3], prhs[4],
              prhs[5], plhs[0], update);
      } else {
        if (is_mixed)
          std::tie(flag, iters, res, tt) =
              hilucsi4m::KSP_solve<true, complex_t>(
                  id, restart, maxit, rtol, verbose, prhs[2], prhs[3], prhs[4],
                  prhs[5], plhs[0], update);
        else
          std::tie(flag, iters, res, tt) =
              hilucsi4m::KSP_solve<false, complex_t>(
                  id, restart, maxit, rtol, verbose, prhs[2], prhs[3], prhs[4],
                  prhs[5], plhs[0], update);
      }
    } else {
      const mxArray* fhdl = nullptr;  // function handler
      const char* feval = "feval";    // feval is used to call user func handler
      char* fname = nullptr;
      int cst_nsp[2];
      const int* nsp_ptr = nullptr;
      if (mxIsClass(prhs[12], "function_handle")) {
        // function handle
        fhdl = prhs[12];
        fname = (char*)feval;
      } else if (mxIsChar(prhs[12])) {
        // function script
        fname = mxArrayToString(prhs[12]);
      } else {
        if (mxGetM(prhs[12]) * mxGetN(prhs[12]) < 2)
          mexErrMsgIdAndTxt("hilucsi4m:mexgateway:ksp_solve",
                            "constant null space mode requires two input "
                            "specifying the range");
        cst_nsp[0] = *mxGetPr(prhs[12]);
        cst_nsp[1] = *(mxGetPr(prhs[12]) + 1);
        nsp_ptr = cst_nsp;
      }
      if (is_real) {
        if (is_mixed)
          std::tie(flag, iters, res, tt) = hilucsi4m::KSP_solve<true>(
              id, restart, maxit, rtol, verbose, prhs[2], prhs[3], prhs[4],
              prhs[5], plhs[0], update, nsp_ptr, fname, fhdl);
        else
          std::tie(flag, iters, res, tt) = hilucsi4m::KSP_solve<false>(
              id, restart, maxit, rtol, verbose, prhs[2], prhs[3], prhs[4],
              prhs[5], plhs[0], update, nsp_ptr, fname, fhdl);
      } else {
        if (is_mixed)
          std::tie(flag, iters, res, tt) =
              hilucsi4m::KSP_solve<true, complex_t>(
                  id, restart, maxit, rtol, verbose, prhs[2], prhs[3], prhs[4],
                  prhs[5], plhs[0], update, nsp_ptr, fname, fhdl);
        else
          std::tie(flag, iters, res, tt) =
              hilucsi4m::KSP_solve<false, complex_t>(
                  id, restart, maxit, rtol, verbose, prhs[2], prhs[3], prhs[4],
                  prhs[5], plhs[0], update, nsp_ptr, fname, fhdl);
      }
      if (fname && fname != feval) mxFree(fname);
    }
    if (nlhs < 2) {
      // handle flag
      if (flag) {
        const char* msg = hif::ksp::flag_repr("FGMRES", flag).c_str();
        mexErrMsgIdAndTxt("hilucsi4m:mexgateway:ksp_solve",
                          "FGMRES failed with flag %d and message:\n%s", flag,
                          msg);
      }
    } else {
      plhs[1] = mxCreateDoubleScalar((double)flag);
      if (nlhs > 2) plhs[2] = mxCreateDoubleScalar((double)iters);
      if (nlhs > 3) plhs[3] = mxCreateDoubleScalar(res);
      if (nlhs > 4) plhs[4] = mxCreateDoubleScalar(tt);
    }
    return;
  }

  if (action == hilucsi4m::HILUCSI4M_KSP_NULL) {
    // KSP null solver
    // act, dbase, rowptr, colind, val, b, (restart, rtol, maxit, x0, verbose,
    // hiprec)
    if (nrhs < 6)
      mexErrMsgIdAndTxt("hilucsi4m:mexgateway:ksp_solve",
                        "KSP solver requires at least 6 inputs");
    if (nlhs < 1)
      mexErrMsgIdAndTxt("hilucsi4m:mexgateway:ksp_solve",
                        "KSP solver requires at least 1 output(s)");
    const int restart = nrhs < 7 ? 30 : (int)mxGetScalar(prhs[6]);
    const double rtol = nrhs < 8 ? 1e-6 : mxGetScalar(prhs[7]);
    const int maxit = nrhs < 9 ? 500 : (int)mxGetScalar(prhs[8]);
    if (nrhs > 9) {
      // has initial guess
      plhs[0] = mxDuplicateArray(prhs[9]);
    } else {
      plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[5]), 1,
                                     is_real ? mxREAL : mxCOMPLEX);
    }
    const bool verbose = nrhs < 11 ? true : (bool)mxGetScalar(prhs[10]);
    const bool hiprec = nrhs < 12 ? false : (bool)mxGetScalar(prhs[11]);
    int flag, iters;
    double res, tt;
    if (is_real) {
      if (is_mixed) {
        if (hiprec)
          std::tie(flag, iters, res, tt) =
              hilucsi4m::KSP_null_solve<true, true>(id, restart, maxit, rtol,
                                                    verbose, prhs[2], prhs[3],
                                                    prhs[4], prhs[5], plhs[0]);
        else
          std::tie(flag, iters, res, tt) =
              hilucsi4m::KSP_null_solve<true, false>(id, restart, maxit, rtol,
                                                     verbose, prhs[2], prhs[3],
                                                     prhs[4], prhs[5], plhs[0]);
      } else {
        if (hiprec)
          std::tie(flag, iters, res, tt) =
              hilucsi4m::KSP_null_solve<false, true>(id, restart, maxit, rtol,
                                                     verbose, prhs[2], prhs[3],
                                                     prhs[4], prhs[5], plhs[0]);
        else
          std::tie(flag, iters, res, tt) =
              hilucsi4m::KSP_null_solve<false, false>(
                  id, restart, maxit, rtol, verbose, prhs[2], prhs[3], prhs[4],
                  prhs[5], plhs[0]);
      }
    } else {
      if (is_mixed) {
        if (hiprec)
          std::tie(flag, iters, res, tt) =
              hilucsi4m::KSP_null_solve<true, true, complex_t>(
                  id, restart, maxit, rtol, verbose, prhs[2], prhs[3], prhs[4],
                  prhs[5], plhs[0]);
        else
          std::tie(flag, iters, res, tt) =
              hilucsi4m::KSP_null_solve<true, false, complex_t>(
                  id, restart, maxit, rtol, verbose, prhs[2], prhs[3], prhs[4],
                  prhs[5], plhs[0]);
      } else {
        if (hiprec)
          std::tie(flag, iters, res, tt) =
              hilucsi4m::KSP_null_solve<false, true, complex_t>(
                  id, restart, maxit, rtol, verbose, prhs[2], prhs[3], prhs[4],
                  prhs[5], plhs[0]);
        else
          std::tie(flag, iters, res, tt) =
              hilucsi4m::KSP_null_solve<false, false, complex_t>(
                  id, restart, maxit, rtol, verbose, prhs[2], prhs[3], prhs[4],
                  prhs[5], plhs[0]);
      }
    }
    if (nlhs < 2) {
      // handle flag
      if (flag != hif::ksp::STAGNATED) {
        const char* msg = hif::ksp::flag_repr("GMRES_Null", flag).c_str();
        mexErrMsgIdAndTxt("hilucsi4m:mexgateway:ksp_solve",
                          "GMRES_Null failed with flag %d and message:\n%s",
                          flag, msg);
      }
    } else {
      plhs[1] = mxCreateDoubleScalar((double)flag);
      if (nlhs > 2) plhs[2] = mxCreateDoubleScalar((double)iters);
      if (nlhs > 3) plhs[3] = mxCreateDoubleScalar(res);
      if (nlhs > 4) plhs[4] = mxCreateDoubleScalar(tt);
    }
    return;
  }

  if (action != hilucsi4m::HILUCSI4M_EXPORT_DATA)
    mexErrMsgIdAndTxt("hilucsi4m:mexgateway", "invalid action flag %d", action);

  // act, dbase, destroy
  const bool destroy = nrhs > 2 ? (bool)mxGetScalar(prhs[2]) : false;
  if (is_real)
    is_mixed ? hilucsi4m::M_export<true>(id, destroy, plhs)
             : hilucsi4m::M_export<false>(id, destroy, plhs);
  else
    is_mixed ? hilucsi4m::M_export<true, complex_t>(id, destroy, plhs)
             : hilucsi4m::M_export<false, complex_t>(id, destroy, plhs);
}

static bool decode_database_struct(const mxArray* rhs, int& id, bool& is_mixed,
                                   bool& is_real) {
  auto id_mx = mxGetField(rhs, 0, fnames[0]);
  if (!id_mx) return true;
  if (!mxIsInt32(id_mx)) return true;
  auto is_mixed_mx = mxGetField(rhs, 0, fnames[1]);
  if (!is_mixed_mx) return true;
  if (!mxIsLogical(is_mixed_mx)) return true;
  auto is_real_mx = mxGetField(rhs, 0, fnames[2]);
  if (!is_real_mx) return true;
  if (!mxIsLogical(is_real_mx)) return true;
  id = mxGetScalar(id_mx);
  is_mixed = *(mxLogical*)mxGetData(is_mixed_mx);
  is_real = *(mxLogical*)mxGetData(is_real_mx);
  return false;
}

template <bool IsMixed, bool isReal>
static void get_fac_info(int id, mwSize& nnz, mwSize& nnz_ef, mwSize& deferrals,
                         mwSize& dy_deferrals, mwSize& drops,
                         mwSize& space_drops, mwSize& lvls) {
  using value_type =
      typename std::conditional<isReal, double, std::complex<double>>::type;
  auto data =
      hilucsi4m::database<IsMixed, value_type>(hilucsi4m::HILUCSI4M_GET, id);
  nnz = data->M->nnz();
  nnz_ef = data->M->nnz_ef();
  deferrals = data->M->stats(0);
  dy_deferrals = data->M->stats(1);
  drops = data->M->stats(4);
  space_drops = data->M->stats(5);
  lvls = data->M->levels();
}
