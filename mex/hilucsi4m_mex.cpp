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
  HILUCSI4M_FACTORIZE = HILUCSI4M_DESTROY + 1,    ///< Factorize
  HILUCSI4M_M_SOLVE   = HILUCSI4M_FACTORIZE + 1,  ///< preconditioner solve
  HILUCSI4M_KSP_SOLVE = HILUCSI4M_M_SOLVE + 1,    ///< KSP solve
};
}

// structure of database fields
static const char* fnames[2] = {"id", "is_mixed"};

// helper function for getting id and is_mixed
static bool decode_database_struct(const mxArray* rhs, int& id, bool& is_mixed);

template <bool IsMixed>
static void get_fac_info(int id, mwSize& nnz, mwSize& nnz_ef, mwSize& deferrals,
                         mwSize& dy_deferrals, mwSize& drops,
                         mwSize& space_drops, mwSize& lvls);

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
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
    int       id;
    if (!is_mixed)
      hilucsi4m::database<false>(hilucsi4m::HILUCSI4M_CREATE, id);
    else
      hilucsi4m::database<true>(hilucsi4m::HILUCSI4M_CREATE, id);
    // create structure
    plhs[0]    = mxCreateStructMatrix(1, 1, 2, fnames);
    auto id_mx = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    *(int*)mxGetData(id_mx) = id;
    mxSetFieldByNumber(plhs[0], 0, 0, id_mx);
    auto is_mixed_mx = mxCreateNumericMatrix(1, 1, mxLOGICAL_CLASS, mxREAL);
    *(mxLogical*)mxGetData(is_mixed_mx) = is_mixed;
    mxSetFieldByNumber(plhs[0], 0, 1, is_mixed_mx);
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
  int  id;
  bool is_mixed;
  if (decode_database_struct(prhs[1], id, is_mixed))
    mexErrMsgIdAndTxt("hilucsi4m:mexgateway",
                      "failed to decode the database structure");
  // destroy
  if (action == hilucsi4m::HILUCSI4M_DESTROY) {
    if (!is_mixed)
      hilucsi4m::database<false>(hilucsi4m::HILUCSI4M_DESTROY, id);
    else
      hilucsi4m::database<true>(hilucsi4m::HILUCSI4M_DESTROY, id);
    return;
  }  // end destroying

  // factorization
  if (action == hilucsi4m::HILUCSI4M_FACTORIZE) {
    // act, dbase, rowptr, colind, val, opt
    if (nrhs < 6)
      mexErrMsgIdAndTxt("hilucsi4m:mexgateway:fac",
                        "factorization requires 6 inputs");
    const double tt = is_mixed ? hilucsi4m::factorize<true>(
                                     id, prhs[2], prhs[3], prhs[4], prhs[5])
                               : hilucsi4m::factorize<false>(
                                     id, prhs[2], prhs[3], prhs[4], prhs[5]);
    if (nlhs) plhs[0] = mxCreateDoubleScalar(tt);
    if (nlhs > 1) {
      static const char* fac_info_names[7] = {
          "levels",       "nnz",   "nnz_EF",     "deferrals",
          "dy_deferrals", "drops", "space_drops"};
      // create an information structure
      mwSize info[7];
      mwSize nnz, nnz_ef, deferrals, dy_deferrals, drops, space_drops, lvls;
      if (!is_mixed)
        get_fac_info<false>(id, info[1], info[2], info[3], info[4], info[5],
                            info[6], info[0]);
      else
        get_fac_info<true>(id, nnz, nnz_ef, deferrals, dy_deferrals, drops,
                           space_drops, lvls);
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
    // act, dbase, b, x
    if (nrhs < 3)
      mexErrMsgIdAndTxt("hilucsi4m:mexgateway:m_solve",
                        "Accessing inv(M) requires 4 inputs");
    // do not use R2017b api for shared copy
    if (nlhs < 1)
      mexErrMsgIdAndTxt("hilucsi4m:mexgateway:m_solve",
                        "Accessing inv(M) requires 1 output");
    plhs[0]         = mxDuplicateArray(prhs[2]);
    const double tt = is_mixed
                          ? hilucsi4m::M_solve<true>(id, prhs[2], plhs[0])
                          : hilucsi4m::M_solve<false>(id, prhs[0], plhs[0]);
    if (nlhs > 1) plhs[1] = mxCreateDoubleScalar(tt);
    return;
  }  // end M solve

  if (action != hilucsi4m::HILUCSI4M_KSP_SOLVE)
    mexErrMsgIdAndTxt("hilucsi4m:mexgateway", "invalid action flag %d", action);

  // KSP solver
  // act, dbase, rowptr, colind, val, b, (restart, rtol, maxit, x0, verbose)
  if (nrhs < 6)
    mexErrMsgIdAndTxt("hilucsi4m:mexgateway:ksp_solve",
                      "KSP solver requires at least 6 inputs");
  if (nlhs < 1)
    mexErrMsgIdAndTxt("hilucsi4m:mexgateway:ksp_solve",
                      "KSP solver requires at least 1 output(s)");
  const int    restart = nrhs < 7 ? 30 : (int)mxGetScalar(prhs[6]);
  const double rtol    = nrhs < 8 ? 1e-6 : mxGetScalar(prhs[7]);
  const int    maxit   = nrhs < 9 ? 500 : (int)mxGetScalar(prhs[8]);
  if (nrhs > 9) {
    // has initial guess
    plhs[0] = mxDuplicateArray(prhs[9]);
  } else {
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[5]), 1, mxREAL);
  }
  const bool verbose = nrhs < 11 ? true : (bool)mxGetScalar(prhs[10]);
  int        flag, iters;
  double     tt;
  if (is_mixed)
    std::tie(flag, iters, tt) =
        hilucsi4m::KSP_solve<true>(id, restart, maxit, rtol, verbose, prhs[2],
                                   prhs[3], prhs[4], prhs[5], plhs[0]);
  else
    std::tie(flag, iters, tt) =
        hilucsi4m::KSP_solve<false>(id, restart, maxit, rtol, verbose, prhs[2],
                                    prhs[3], prhs[4], prhs[5], plhs[0]);
  if (nlhs < 2) {
    // handle flag
    if (flag) {
      const char* msg = hilucsi::ksp::flag_repr("FGMRES", flag).c_str();
      mexErrMsgIdAndTxt("hilucsi4m:mexgateway:ksp_solve",
                        "FGMRES failed with flag %d and message:\n%s", flag,
                        msg);
    }
  } else {
    plhs[1] = mxCreateDoubleScalar((double)flag);
    if (nlhs > 2) plhs[2] = mxCreateDoubleScalar((double)iters);
    if (nlhs > 3) plhs[3] = mxCreateDoubleScalar(tt);
  }
}

static bool decode_database_struct(const mxArray* rhs, int& id,
                                   bool& is_mixed) {
  auto id_mx = mxGetField(rhs, 0, fnames[0]);
  if (!id_mx) return true;
  if (!mxIsInt32(id_mx)) return true;
  auto is_mixed_mx = mxGetField(rhs, 0, fnames[1]);
  if (!is_mixed_mx) return true;
  if (!mxIsLogical(is_mixed_mx)) return true;
  id       = mxGetScalar(id_mx);
  is_mixed = *(mxLogical*)mxGetData(is_mixed_mx);
  return false;
}

template <bool IsMixed>
static void get_fac_info(int id, mwSize& nnz, mwSize& nnz_ef, mwSize& deferrals,
                         mwSize& dy_deferrals, mwSize& drops,
                         mwSize& space_drops, mwSize& lvls) {
  auto data    = hilucsi4m::database<IsMixed>(hilucsi4m::HILUCSI4M_GET, id);
  nnz          = data->M->nnz();
  nnz_ef       = data->M->nnz_EF();
  deferrals    = data->M->stats(0);
  dy_deferrals = data->M->stats(1);
  drops        = data->M->stats(4);
  space_drops  = data->M->stats(5);
  lvls         = data->M->levels();
}
