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

#include <complex>
#include <type_traits>

#include "mex.h"

#ifdef HIFIR4M_USE_32INT
typedef int integer_type;
#define INT_TYPE mxINT32_CLASS
#else
typedef mwSignedIndex integer_type;
#define INT_TYPE mxINT64_CLASS
#endif

template <bool IsReal>
static inline typename std::enable_if<IsReal, void>::type convert_ijv2crs(
    const mwSize n, const mwSize nnz, const integer_type *rows,
    const integer_type *cols, const mxArray *vs, integer_type *rowptr,
    integer_type *colind, mxArray *vals) {
  // handle Fortran indices
  for (mwSize i(0); i < nnz; ++i)
    ++rowptr[rows[i]];  // implicit +1 per Fortran index
  rowptr[0] = 0;
  for (mwSize i(0); i < n; ++i) rowptr[i + 1] += rowptr[i];
  --rowptr;
  double *vals_ptr = mxGetPr(vals), *vs_ptr = (double *)mxGetPr(vs);
  for (mwSize i(0); i < nnz; ++i) {
    const auto idx = rowptr[rows[i]];
    vals_ptr[idx] = vs_ptr[i];
    colind[idx] = cols[i] - 1;
    ++rowptr[rows[i]];
  }
  // revert rowptr (notice that the address has been shiftted leftward by 1)
  for (mwSize i(n); i > 1u; --i) rowptr[i] = rowptr[i - 1];
  rowptr[1] = 0;
}

template <bool IsReal>
static inline typename std::enable_if<!IsReal, void>::type convert_ijv2crs(
    const mwSize n, const mwSize nnz, const integer_type *rows,
    const integer_type *cols, const mxArray *vs, integer_type *rowptr,
    integer_type *colind, mxArray *vals) {
  // complex arithmetic

  // handle Fortran indices
  for (mwSize i(0); i < nnz; ++i)
    ++rowptr[rows[i]];  // implicit +1 per Fortran index
  rowptr[0] = 0;
  for (mwSize i(0); i < n; ++i) rowptr[i + 1] += rowptr[i];
  --rowptr;
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && MX_HAS_INTERLEAVED_COMPLEX
  std::complex<double> *vals_ptr = (std::complex<double> *)mxGetData(vals),
                       *vs_ptr = (std::complex<double> *)mxGetData(vs);
#else
  double *vals_r_ptr = mxGetPr(vals), *vals_i_ptr = mxGetPi(vals),
         *vs_r_ptr = (double *)mxGetPr(vs), *vs_i_ptr = mxGetPi(vs);
#endif
  for (mwSize i(0); i < nnz; ++i) {
    const auto idx = rowptr[rows[i]];
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && MX_HAS_INTERLEAVED_COMPLEX
    vals_ptr[idx] = vs_ptr[i];
#else
    vals_r_ptr[idx] = vs_r_ptr[i];
    vals_i_ptr[idx] = vs_i_ptr[i];
#endif
    colind[idx] = cols[i] - 1;
    ++rowptr[rows[i]];
  }
  // revert rowptr (notice that the address has been shiftted leftward by 1)
  for (mwSize i(n); i > 1u; --i) rowptr[i] = rowptr[i - 1];
  rowptr[1] = 0;
}

// convert IJV format from MATLAB built-in find to CRS
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // n, I, J, V
  if (nrhs < 4)
    mexErrMsgIdAndTxt("hilucsi4m:ijv2crs", "this routine requires 4 inputs");

  if (nlhs != 3)
    mexErrMsgIdAndTxt("hilucsi4m:ijv2crs", "this routine requires 3 outputs");

  const mwSize n = mxGetScalar(prhs[0]);
  const mwSize nnz = mxGetM(prhs[1]);
  if (nnz != mxGetM(prhs[2]) || nnz != mxGetM(prhs[3]))
    mexErrMsgIdAndTxt("hilucsi4m:ijv2crs", "inconsistent NNZ for IJV");

  plhs[0] = mxCreateNumericMatrix(n + 1, 1, INT_TYPE, mxREAL);
  plhs[1] = mxCreateUninitNumericMatrix(nnz, 1, INT_TYPE, mxREAL);
  if (mxIsComplex(prhs[3]))
    plhs[2] = mxCreateDoubleMatrix(nnz, 1, mxCOMPLEX);
  else
    plhs[2] = mxCreateDoubleMatrix(nnz, 1, mxREAL);

  // convert
  if (!mxIsComplex(prhs[3]))
    convert_ijv2crs<true>(n, nnz, (const integer_type *)mxGetData(prhs[1]),
                          (const integer_type *)mxGetData(prhs[2]), prhs[3],
                          (integer_type *)mxGetData(plhs[0]),
                          (integer_type *)mxGetData(plhs[1]), plhs[2]);
  else
    convert_ijv2crs<false>(n, nnz, (const integer_type *)mxGetData(prhs[1]),
                           (const integer_type *)mxGetData(prhs[2]), prhs[3],
                           (integer_type *)mxGetData(plhs[0]),
                           (integer_type *)mxGetData(plhs[1]), plhs[2]);
}
