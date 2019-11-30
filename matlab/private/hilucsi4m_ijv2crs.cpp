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

#include "mex.h"

static void convert_ijv2crs(const mwSize n, const mwSize nnz, const int *rows,
                            const int *cols, const double *vs, int *rowptr,
                            int *colind, double *vals);

// convert IJV format from MATLAB built-in find to CRS
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // n, I, J, V
  if (nrhs < 4)
    mexErrMsgIdAndTxt("hilucsi4m:ijv2crs", "this routine requires 4 inputs");

  if (nlhs != 3)
    mexErrMsgIdAndTxt("hilucsi4m:ijv2crs", "this routine requires 3 outputs");

  const mwSize n   = mxGetScalar(prhs[0]);
  const mwSize nnz = mxGetM(prhs[1]);
  if (nnz != mxGetM(prhs[2]) || nnz != mxGetM(prhs[3]))
    mexErrMsgIdAndTxt("hilucsi4m:ijv2crs", "inconsistent NNZ for IJV");

  plhs[0] = mxCreateNumericMatrix(n + 1, 1, mxINT32_CLASS, mxREAL);
  plhs[1] = mxCreateUninitNumericMatrix(nnz, 1, mxINT32_CLASS, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(nnz, 1, mxREAL);

  // convert
  convert_ijv2crs(n, nnz, (const int *)mxGetData(prhs[1]),
                  (const int *)mxGetData(prhs[2]), mxGetPr(prhs[3]),
                  (int *)mxGetData(plhs[0]), (int *)mxGetData(plhs[1]),
                  mxGetPr(plhs[2]));
}

static void convert_ijv2crs(const mwSize n, const mwSize nnz, const int *rows,
                            const int *cols, const double *vs, int *rowptr,
                            int *colind, double *vals) {
  // handle Fortran indices
  for (mwSize i(0); i < nnz; ++i)
    ++rowptr[rows[i]];  // implicit +1 per Fortran index
  rowptr[0] = 0;
  for (mwSize i(0); i < n; ++i) rowptr[i + 1] += rowptr[i];
  --rowptr;
  for (mwSize i(0); i < nnz; ++i) {
    const auto idx = rowptr[rows[i]];
    vals[idx]      = vs[i];
    colind[idx]    = cols[i] - 1;
    ++rowptr[rows[i]];
  }
  // revert rowptr (notice that the address has been shiftted leftward by 1)
  for (mwSize i(n); i > 1u; --i) rowptr[i] = rowptr[i - 1];
  rowptr[1] = 0;
}
