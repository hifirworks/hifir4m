/*
                This file is part of HILUCSI4M project

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

#include <complex>

#include "mex.h"

// convert IJV format from MATLAB built-in find to CRS
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  (void)nlhs;
  (void)nrhs;
  (void)prhs;
  const bool is_int64 =
#ifdef HIFIR4M_USE_32INT
      false;
#else
      true;
#endif
  plhs[0] = mxCreateLogicalScalar(is_int64);
}
