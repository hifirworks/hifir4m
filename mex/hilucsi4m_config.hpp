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

// WARNING! This file must be included before including HILUCSI.hpp

#ifndef HILUCSI4M_CONFIG_HPP_
#define HILUCSI4M_CONFIG_HPP_

// mex
#include "mex.h"

// configure the printing to matlab command window stead of stdout/stderr
#ifdef HILUCSI_STDOUT
#  undef HILUCSI_STDOUT
#endif
#ifdef HILUCSI_STDERR
#  undef HILUCSI_STDERR
#endif
#define HILUCSI_STDOUT(__msg) \
  do {                        \
    mexPrintf(__msg);         \
  } while (false)
#define HILUCSI_STDERR(__msg) \
  do {                        \
    mexWarnMsgTxt(__msg);     \
  } while (false)

// enforce throw instead of calling abort!
#ifndef HILUCSI_THROW
#  define HILUCSI_THROW
#endif  // HILUCSI_THROW

// lapack integer
#ifdef HILUCSI_LAPACK_INT
#  undef HILUCSI_LAPACK_INT
#endif
#define HILUCSI_LAPACK_INT mwSignedIndex

#endif  // HILUCSI4M_CONFIG_HPP_
