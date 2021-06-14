/*
                This file is part of HIFIR4M project

    Copyright (C) 2021 NumGeom Group at Stony Brook University

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

#ifndef HIFIR4M_NUMARRAY_HPP_
#define HIFIR4M_NUMARRAY_HPP_

#include <complex>
#include <vector>

#include "mex.h"

namespace hifir4m {

/*!
 * @brief Get the complex data from a complex mex array
 * @param[in] v Input mex array
 * @param[out] ptr Output pointer to the proper complex data
 * @param[out] buf Buffer space if no interleaved complex is supported
 */
template <class T>
inline void get_num_data(const mxArray *v, std::complex<T> **ptr,
                         std::vector<std::complex<T>> &buf) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && MX_HAS_INTERLEAVED_COMPLEX
  // if we have interleaved complex arrays, then directly retrieving the pointer
  (void)buf;
  *ptr = (std::complex<T> *)mxGetData(v);
#else
  // if we don't have interleaved complex arrays, use buffer
  const auto n = mxGetM(v);
  buf.resize(n);
  T *pr = (T *)mxGetData(v);
  T *pi = (T *)mxGetImagData(v);
  for (std::size_t i(0); i < n; ++i) {
    buf[i].real(pr[i]);
    buf[i].imag(pi[i]);
  }
  *ptr = buf.data();
#endif
}

// overload for double
template <class T>
inline void get_num_data(const mxArray *v, T **ptr, std::vector<T> &) {
  *ptr = (T *)mxGetData(v);
}

/*!
 * @brief Set the complex data to mex array
 */
template <class T>
inline void set_num_data(mxArray *v, const std::complex<T> *ptr) {
#if defined(MX_HAS_INTERLEAVED_COMPLEX) && MX_HAS_INTERLEAVED_COMPLEX
  // NOTE if we have interleaved complex arrays, then ptr should point
  // to v, thus nothing needs to be done.
  (void)v;
  (void)ptr;
#else
  const auto n = mxGetM(v);
  T *pr = (T *)mxGetData(v);
  T *pi = (T *)mxGetImagData(v);
  for (std::size_t i(0); i < n; ++i) {
    pr[i] = ptr[i].real();
    pi[i] = ptr[i].imag();
  }
#endif
}

// overload for doubles
template <class T>
inline void set_num_data(mxArray *, const T *) {}

}  // namespace hifir4m

#endif  // HIFIR4M_NUMARRAY_HPP_