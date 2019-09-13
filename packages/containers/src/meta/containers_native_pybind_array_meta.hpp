/*
//@HEADER
// ************************************************************************
//
// containers_native_pybind_array_meta.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifdef HAVE_PYBIND11
#ifndef CONTAINERS_NATIVE_PYBIND11_ARRAY_HPP_
#define CONTAINERS_NATIVE_PYBIND11_ARRAY_HPP_

#include "containers_meta_basic.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>


namespace pressio{ namespace containers{ namespace meta {

/*
 * this metafunction is here because a pybind11::array_t
 * can have arbitrary size sine it maps to numpy.
 * And I don't know yet if a pybind11::array_t can be
 * checked at compile time to be a vector or matrix
 */

template <typename T, typename enable = void>
struct is_cstyle_array_pybind11 : std::false_type {};

template <typename T>
struct is_cstyle_array_pybind11<
  T,
  ::pressio::mpl::enable_if_t<
    mpl::is_same<
      T,
      pybind11::array_t<
	typename T::value_type,
	pybind11::array::c_style
	>
      >::value
    >
  > : std::true_type{};
//----------------------------------------------

template <typename T, typename enable = void>
struct is_array_pybind11 : std::false_type {};

template <typename T>
struct is_array_pybind11<
  T,
  ::pressio::mpl::enable_if_t<
    is_cstyle_array_pybind11<T>::value
    >
  > : std::true_type{};


}}}//end namespace pressio::containers::meta
#endif
#endif
