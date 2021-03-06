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

#include "hifir4m.hpp"

// structure of database fields
static const char* fnames[3] = {"id", "is_mixed", "is_real"};

// structure of HIF information
static const char* fac_info_names[11] = {
    "nnz",         "nnz_ef", "nnz_ldu", "deferrals",  "dy_deferrals", "drops",
    "space_drops", "levels", "rank",    "schur_rank", "schur_size"};

// helper function for getting id and is_mixed
static bool decode_database_struct(const mxArray* rhs, int& id, bool& is_mixed,
                                   bool& is_real);

template <bool IsMixed, bool isReal>
static void get_fac_info(int id, mwSize& nnz, mwSize& nnz_ef, mwSize& nnz_ldu,
                         mwSize& deferrals, mwSize& dy_deferrals, mwSize& drops,
                         mwSize& space_drops, mwSize& lvls, mwSize& rank,
                         mwSize& schur_rank, mwSize& schur_size);

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  using complex_t = std::complex<double>;

  if (nrhs < 1)
    mexErrMsgIdAndTxt("hifir4m:mexgateway",
                      "At least one argument must be given");
  const int action = (int)mxGetScalar(prhs[0]);
  // create
  if (action == hifir4m::HIFIR4M_CREATE) {
    // determine if use mixed precision
    const int is_mixed = nrhs < 2 ? 0 : (int)mxGetScalar(prhs[1]);
    const int is_complex = nrhs < 3 ? 0 : (int)mxGetScalar(prhs[2]);
    int id;
    if (!is_mixed) {
      if (!is_complex)
        hifir4m::database<false, double>(hifir4m::HIFIR4M_CREATE, id);
      else
        hifir4m::database<false, complex_t>(hifir4m::HIFIR4M_CREATE, id);
    } else {
      if (!is_complex)
        hifir4m::database<true, double>(hifir4m::HIFIR4M_CREATE, id);
      else
        hifir4m::database<true, complex_t>(hifir4m::HIFIR4M_CREATE, id);
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
    mexErrMsgIdAndTxt("hifir4m:mexgateway",
                      "the database token structure must be passed in");
  if (!mxIsStruct(prhs[1]))
    mexErrMsgIdAndTxt("hifir4m:mexgateway",
                      "the database token must be structure");
  int id;
  bool is_mixed, is_real;
  if (decode_database_struct(prhs[1], id, is_mixed, is_real))
    mexErrMsgIdAndTxt("hifir4m:mexgateway",
                      "failed to decode the database structure");
  // clear database
  if (action == hifir4m::HIFIR4M_CLEAR) {
    if (!is_mixed) {
      if (is_real)
        hifir4m::database<false, double>(hifir4m::HIFIR4M_CLEAR, id);
      else
        hifir4m::database<false, complex_t>(hifir4m::HIFIR4M_CLEAR, id);
    } else {
      if (is_real)
        hifir4m::database<true, double>(hifir4m::HIFIR4M_CLEAR, id);
      else
        hifir4m::database<true, complex_t>(hifir4m::HIFIR4M_CLEAR, id);
    }
    return;
  }  // end clearing

  if (action == hifir4m::HIFIR4M_CHECK) {
    if (!is_mixed) {
      if (is_real) {
        auto* dbase = hifir4m::database<false>(hifir4m::HIFIR4M_GET, id);
        if (!dbase || !dbase->M)
          plhs[0] = mxCreateLogicalScalar(1);
        else
          plhs[0] = mxCreateLogicalScalar((mxLogical)dbase->M->empty());
      } else {
        auto* dbase =
            hifir4m::database<false, complex_t>(hifir4m::HIFIR4M_GET, id);
        if (!dbase || !dbase->M)
          plhs[0] = mxCreateLogicalScalar(1);
        else
          plhs[0] = mxCreateLogicalScalar((mxLogical)dbase->M->empty());
      }
    } else {
      if (is_real) {
        auto* dbase = hifir4m::database<true>(hifir4m::HIFIR4M_GET, id);
        if (!dbase || !dbase->M)
          plhs[0] = mxCreateLogicalScalar(1);
        else
          plhs[0] = mxCreateLogicalScalar((mxLogical)dbase->M->empty());
      } else {
        auto* dbase =
            hifir4m::database<true, complex_t>(hifir4m::HIFIR4M_GET, id);
        if (!dbase || !dbase->M)
          plhs[0] = mxCreateLogicalScalar(1);
        else
          plhs[0] = mxCreateLogicalScalar((mxLogical)dbase->M->empty());
      }
    }
    return;
  }

  // destroy
  if (action == hifir4m::HIFIR4M_DESTROY) {
    if (!is_mixed) {
      if (is_real)
        hifir4m::database<false>(hifir4m::HIFIR4M_DESTROY, id);
      else
        hifir4m::database<false, complex_t>(hifir4m::HIFIR4M_DESTROY, id);
    } else {
      if (is_real)
        hifir4m::database<true>(hifir4m::HIFIR4M_DESTROY, id);
      else
        hifir4m::database<true, complex_t>(hifir4m::HIFIR4M_DESTROY, id);
    }
    return;
  }  // end destroying

  // factorization
  if (action == hifir4m::HIFIR4M_FACTORIZE) {
    // act, dbase, rowptr, colind, val, opt
    if (nrhs < 6)
      mexErrMsgIdAndTxt("hifir4m:mexgateway:fac",
                        "factorization requires 6 inputs");
    double tt;
    if (is_real)
      tt = is_mixed ? hifir4m::factorize<true>(id, prhs[2], prhs[3], prhs[4],
                                               prhs[5])
                    : hifir4m::factorize<false>(id, prhs[2], prhs[3], prhs[4],
                                                prhs[5]);
    else
      tt = is_mixed ? hifir4m::factorize<true, complex_t>(id, prhs[2], prhs[3],
                                                          prhs[4], prhs[5])
                    : hifir4m::factorize<false, complex_t>(id, prhs[2], prhs[3],
                                                           prhs[4], prhs[5]);
    // create an information structure
    mwSize info[11];
    if (!is_mixed) {
      if (is_real)
        get_fac_info<false, true>(id, info[0], info[1], info[2], info[3],
                                  info[4], info[5], info[6], info[7], info[8],
                                  info[9], info[10]);
      else
        get_fac_info<false, false>(id, info[0], info[1], info[2], info[3],
                                   info[4], info[5], info[6], info[7], info[8],
                                   info[9], info[10]);
    } else {
      if (is_real)
        get_fac_info<true, true>(id, info[0], info[1], info[2], info[3],
                                 info[4], info[5], info[6], info[7], info[8],
                                 info[9], info[10]);
      else
        get_fac_info<true, false>(id, info[0], info[1], info[2], info[3],
                                  info[4], info[5], info[6], info[7], info[8],
                                  info[9], info[10]);
    }
    plhs[0] = mxCreateStructMatrix(1, 1, 11, fac_info_names);
    for (int i = 0; i < 11; ++i) {
      auto mx_ptr = mxCreateUninitNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
      *(std::size_t*)mxGetData(mx_ptr) = info[i];
      mxSetFieldByNumber(plhs[0], 0, i, mx_ptr);
    }
    if (nlhs > 1) plhs[1] = mxCreateDoubleScalar(tt);
    return;
  }  // end factorization

  // query information
  if (action == hifir4m::HIFIR4M_QUERY) {
    mwSize info[11];
    if (!is_mixed) {
      if (is_real)
        get_fac_info<false, true>(id, info[0], info[1], info[2], info[3],
                                  info[4], info[5], info[6], info[7], info[8],
                                  info[9], info[10]);
      else
        get_fac_info<false, false>(id, info[0], info[1], info[2], info[3],
                                   info[4], info[5], info[6], info[7], info[8],
                                   info[9], info[10]);
    } else {
      if (is_real)
        get_fac_info<true, true>(id, info[0], info[1], info[2], info[3],
                                 info[4], info[5], info[6], info[7], info[8],
                                 info[9], info[10]);
      else
        get_fac_info<true, false>(id, info[0], info[1], info[2], info[3],
                                  info[4], info[5], info[6], info[7], info[8],
                                  info[9], info[10]);
    }
    plhs[0] = mxCreateStructMatrix(1, 1, 11, fac_info_names);
    for (int i = 0; i < 11; ++i) {
      auto mx_ptr = mxCreateUninitNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
      *(std::size_t*)mxGetData(mx_ptr) = info[i];
      mxSetFieldByNumber(plhs[0], 0, i, mx_ptr);
    }
    return;
  }

  // M solve
  if (action == hifir4m::HIFIR4M_M_SOLVE) {
    // act, dbase, b, (trans), (r)
    if (nrhs < 3)
      mexErrMsgIdAndTxt("hifir4m:mexgateway:m_solve",
                        "Accessing inv(M) requires 3 inputs");
    double tt;
    if (nrhs <= 5) {
      plhs[0] = mxDuplicateArray(prhs[2]);
      const bool trans = nrhs <= 3 ? false : (bool)mxGetScalar(prhs[3]);
      const std::size_t r = nrhs <= 4 ? 0u : (std::size_t)mxGetScalar(prhs[4]);
      if (is_real)
        tt = is_mixed ? hifir4m::M_solve<true>(id, prhs[2], plhs[0], trans, r)
                      : hifir4m::M_solve<false>(id, prhs[2], plhs[0], trans, r);
      else
        tt = is_mixed ? hifir4m::M_solve<true, complex_t>(id, prhs[2], plhs[0],
                                                          trans, r)
                      : hifir4m::M_solve<false, complex_t>(id, prhs[2], plhs[0],
                                                           trans, r);
    } else {
      // inner iterations
      // act, dbase, b, rowptr, colind, val, N, (trans), (r)
      if (nrhs < 7)
        mexErrMsgIdAndTxt("hifir4m:mexgateway:m_solve_inner",
                          "Accessing inv(M) requires at least 7 inputs");
      const int N = (int)mxGetScalar(prhs[6]);
      const bool trans = nrhs <= 7 ? false : (bool)mxGetScalar(prhs[7]);
      const std::size_t r =
          nrhs <= 8 ? (std::size_t)-1 : (std::size_t)mxGetScalar(prhs[8]);
      plhs[0] = mxDuplicateArray(prhs[2]);
      if (is_real)
        tt = is_mixed ? hifir4m::M_solve<true>(id, prhs[3], prhs[4], prhs[5],
                                               prhs[2], N, plhs[0], trans, r)
                      : hifir4m::M_solve<false>(id, prhs[3], prhs[4], prhs[5],
                                                prhs[2], N, plhs[0], trans, r);
      else
        tt = is_mixed ? hifir4m::M_solve<true, complex_t>(id, prhs[3], prhs[4],
                                                          prhs[5], prhs[2], N,
                                                          plhs[0], trans, r)
                      : hifir4m::M_solve<false, complex_t>(id, prhs[3], prhs[4],
                                                           prhs[5], prhs[2], N,
                                                           plhs[0], trans, r);
    }
    if (nlhs > 1) plhs[1] = mxCreateDoubleScalar(tt);
    return;
  }  // end M solve

  // M solve
  if (action == hifir4m::HIFIR4M_M_MULTIPLY) {
    // act, dbase, b, (trans), (r)
    if (nrhs < 3)
      mexErrMsgIdAndTxt("hifir4m:mexgateway:m_multiply",
                        "M*x requires 3 inputs");
    double tt;
    plhs[0] = mxDuplicateArray(prhs[2]);
    const bool trans = nrhs <= 3 ? false : (bool)mxGetScalar(prhs[3]);
    const std::size_t r = nrhs <= 4 ? 0u : (std::size_t)mxGetScalar(prhs[4]);
    if (is_real)
      tt = is_mixed
               ? hifir4m::M_multiply<true>(id, prhs[2], plhs[0], trans, r)
               : hifir4m::M_multiply<false>(id, prhs[2], plhs[0], trans, r);
    else
      tt = is_mixed ? hifir4m::M_multiply<true, complex_t>(id, prhs[2], plhs[0],
                                                           trans, r)
                    : hifir4m::M_multiply<false, complex_t>(id, prhs[2],
                                                            plhs[0], trans, r);
    if (nlhs > 1) plhs[1] = mxCreateDoubleScalar(tt);
    return;
  }  // end multiply

  if (action != hifir4m::HIFIR4M_EXPORT_DATA)
    mexErrMsgIdAndTxt("hifir4m:mexgateway", "invalid action flag %d", action);

  // act, dbase, destroy
  const bool destroy = nrhs > 2 ? (bool)mxGetScalar(prhs[2]) : false;
  if (is_real)
    is_mixed ? hifir4m::M_export<true>(id, destroy, plhs)
             : hifir4m::M_export<false>(id, destroy, plhs);
  else
    is_mixed ? hifir4m::M_export<true, complex_t>(id, destroy, plhs)
             : hifir4m::M_export<false, complex_t>(id, destroy, plhs);
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
static void get_fac_info(int id, mwSize& nnz, mwSize& nnz_ef, mwSize& nnz_ldu,
                         mwSize& deferrals, mwSize& dy_deferrals, mwSize& drops,
                         mwSize& space_drops, mwSize& lvls, mwSize& rank,
                         mwSize& schur_rank, mwSize& schur_size) {
  using value_type =
      typename std::conditional<isReal, double, std::complex<double>>::type;
  auto data = hifir4m::database<IsMixed, value_type>(hifir4m::HIFIR4M_GET, id);
  if (data->M->empty()) {
    nnz = nnz_ef = nnz_ldu = deferrals = dy_deferrals = drops = space_drops =
        lvls = rank = schur_rank = schur_size = 0u;
    return;
  }
  nnz = data->M->nnz();
  nnz_ef = data->M->nnz_ef();
  nnz_ldu = data->M->nnz_ldu();
  deferrals = data->M->stats(0);
  dy_deferrals = data->M->stats(1);
  drops = data->M->stats(4);
  space_drops = data->M->stats(5);
  lvls = data->M->levels();
  rank = data->M->rank();
  schur_rank = data->M->schur_rank();
  schur_size = data->M->schur_size();
}
