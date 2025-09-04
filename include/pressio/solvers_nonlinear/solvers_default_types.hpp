/*
//@HEADER
// ************************************************************************
//
// solvers_default_types.hpp
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

#ifndef PRESSIOROM_SOLVERS_NONLINEAR_SOLVERS_DEFAULT_TYPES_HPP_
#define PRESSIOROM_SOLVERS_NONLINEAR_SOLVERS_DEFAULT_TYPES_HPP_

namespace pressio{
namespace nlsol{

template<class T, class = void>
struct normal_eqs_default_types{
  using hessian_type  = void;
  using gradient_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct normal_eqs_default_types<
  T, std::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using hessian_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
  using gradient_type = T;

  static hessian_type createHessian(const T & v){
    const auto ext = ::pressio::ops::extent(v, 0);
    return hessian_type(ext, ext);
  }
};
#endif

template<class T> using normal_eqs_default_hessian_t =
  typename normal_eqs_default_types<T>::hessian_type;
template<class T> using normal_eqs_default_gradient_t =
  typename normal_eqs_default_types<T>::gradient_type;

// ====================================================================

template<class T, class = void>
struct valid_state_for_least_squares : std::false_type{};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct valid_state_for_least_squares<
  T, std::enable_if_t< ::pressio::is_vector_eigen<T>::value >
  > : std::true_type{};
#endif

}}
#endif  // PRESSIOROM_SOLVERS_NONLINEAR_SOLVERS_DEFAULT_TYPES_HPP_
