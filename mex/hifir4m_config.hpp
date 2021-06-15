/*
                This file is part of HIFIR4M project

    Copyright (C) 2019--2021 NumGeom Group at Stony Brook University

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

// WARNING! This file must be included before including HIFIR.hpp

#ifndef HIFIR4M_CONFIG_HPP_
#define HIFIR4M_CONFIG_HPP_

// mex
#include "mex.h"

// configure the printing to matlab command window stead of stdout/stderr
#ifdef HIF_STDOUT
#undef HIF_STDOUT
#endif
#ifdef HIF_STDERR
#undef HIF_STDERR
#endif

#if !defined(HAVE_OCTAVE) || !defined(__unix__)
#define ioFlush() mexEvalString("drawnow;")
#else
#define ioFlush()
#endif
#define HIF_STDOUT(__msg) \
  do {                    \
    mexPrintf(__msg);     \
    mexPrintf("\n");      \
    ioFlush(); \
  } while (false)
#define HIF_STDERR(__msg) \
  do {                    \
    mexWarnMsgTxt(__msg); \
    mexWarnMsgTxt("\n");  \
  } while (false)

// No ASCII color
#ifndef HIF_LOG_PLAIN_PREFIX
#define HIF_LOG_PLAIN_PREFIX
#endif

// enforce throw instead of calling abort!
#ifndef HIF_THROW
#define HIF_THROW
#endif  // HIF_THROW

// lapack integer
#ifdef HIF_LAPACK_INT
#undef HIF_LAPACK_INT
#endif

#ifdef HAVE_OCTAVE
#define HIF_LAPACK_INT int
#else
#define HIF_LAPACK_INT mwSignedIndex
#endif

#endif  // HIFIR4M_CONFIG_HPP_
