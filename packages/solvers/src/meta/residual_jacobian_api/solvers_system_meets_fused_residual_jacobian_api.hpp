/*
//@HEADER
// ************************************************************************
//
// solvers_system_meets_fused_residual_jacobian_api.hpp
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

#ifndef SOLVERS_SYSTEM_MEETS_FUSED_RESIDUAL_JACOBIAN_API_HPP_
#define SOLVERS_SYSTEM_MEETS_FUSED_RESIDUAL_JACOBIAN_API_HPP_

namespace pressio{ namespace solvers{ namespace meta {

template<typename T, typename enable = void>
struct system_meets_fused_residual_jacobian_api : std::false_type{};

template<typename T>
struct system_meets_fused_residual_jacobian_api
<T,
 ::pressio::mpl::enable_if_t<
   ::pressio::mpl::is_detected<has_scalar_typedef, T>::value   and
   ::pressio::mpl::is_detected<has_state_typedef, T>::value    and
   ::pressio::mpl::is_detected<has_residual_typedef, T>::value and
   ::pressio::mpl::is_detected<has_jacobian_typedef, T>::value and

   // --- detect createResidualObject ---
   ::pressio::mpl::is_same<
     typename T::residual_type,
     decltype(
	      std::declval<T const>().createResidualObject
	      (
	       std::declval<typename T::state_type const&>()
	       )
	      )
     >::value and

   // --- detect createJacobianObject ---
   ::pressio::mpl::is_same<
     typename T::jacobian_type,
     decltype(
	      std::declval<T const>().createJacobianObject
	      (
	       std::declval<typename T::state_type const&>()
	       )
	      )
     >::value and

   // --- detect residualAndJacobian ---
   std::is_void<
     decltype(
	      std::declval<T const>().residualAndJacobian
	      (
	       std::declval<typename T::state_type const&>(),
	       std::declval<typename T::residual_type &>(),
	       std::declval<typename T::jacobian_type &>(),
	       // the norm type and norm value
	       ::pressio::solvers::Norm::Undefined,
	       std::declval<typename T::scalar_type &>()
	       )
	      )
     >::value and

   // --- detect residualNorm ---
   std::is_void<
     decltype(
	      std::declval<T const>().residualNorm
	      (
	       std::declval<typename T::state_type const &>(),
	       // the norm type and norm value
	       ::pressio::solvers::Norm::Undefined,
	       std::declval<typename T::scalar_type &>()
	       )
	      )
     >::value

   >
 > : std::true_type{};


// template<typename T>
// struct system_meets_fused_residual_jacobian_api
// <T,
//  ::pressio::mpl::enable_if_t<
//    ::pressio::mpl::is_detected<has_scalar_typedef, T>::value   and
//    ::pressio::mpl::is_detected<has_state_typedef, T>::value    and
//    ::pressio::mpl::is_detected<has_residual_typedef, T>::value and
//    ::pressio::mpl::is_detected<has_jacobian_typedef, T>::value and

//    // --- detect residual method with one arg ---
//    ::pressio::mpl::is_same<
//      typename T::residual_type,
//      decltype(
// 	      std::declval<T const>().residual
// 	      ( std::declval<typename T::state_type const&>() )
// 	      )
//      >::value and

//    // --- detect jacobian method with one arg ---
//    ::pressio::mpl::is_same<
//      typename T::jacobian_type,
//      decltype(
// 	      std::declval<T const>().jacobian
// 	      (std::declval<typename T::state_type const&>())
// 	      )
//      >::value and

//    // --- detect residual and jacobian ---
//    std::is_void<
//      decltype(
// 	      std::declval<T const>().residualAndJacobian
// 	      (
// 	       std::declval<typename T::state_type const&>(),
// 	       std::declval<typename T::residual_type &>(),
// 	       std::declval<typename T::jacobian_type &>()
// 	       )
// 	      )
//      >::value
//    >
//  > : std::true_type{};

}}} // namespace pressio::solvers::meta
#endif
