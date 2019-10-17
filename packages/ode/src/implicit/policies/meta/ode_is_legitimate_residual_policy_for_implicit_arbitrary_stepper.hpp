/*
//@HEADER
// ************************************************************************
//
// ode_is_legitimate_residual_policy_for_implicit_arbitrary_stepper.hpp
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

#ifndef ODE_IS_LEGITIMATE_RESIDUAL_POLICY_FOR_IMPLICIT_ARBITRARY_STEPPER_HPP_
#define ODE_IS_LEGITIMATE_RESIDUAL_POLICY_FOR_IMPLICIT_ARBITRARY_STEPPER_HPP_

#include "../base/ode_implicit_residual_policy_base.hpp"
#include "../../../meta/ode_is_legitimate_implicit_state_type.hpp"
#include "../../../meta/ode_is_legitimate_implicit_residual_type.hpp"
#include "../../../meta/ode_has_num_aux_states_static_member.hpp"
#include "../../../meta/ode_has_stepper_order_static_member.hpp"

namespace pressio{ namespace ode{ namespace meta {

template<
  typename T, typename state_t, typename residual_t, typename system_t, typename scalar_t,
  typename enable = void
  >
struct residual_policy_callable_with_six_args : std::false_type{};

template<
  typename T, typename state_t, typename residual_t, typename system_t, typename scalar_t
  >
struct residual_policy_callable_with_six_args<
  T, state_t, residual_t, system_t, scalar_t,
  ::pressio::mpl::enable_if_t<
    std::is_same<
      residual_t,
      decltype
      (
       std::declval<T>().operator()
       ( std::declval<const state_t &>(),
	 std::declval<const std::array<state_t, T::num_aux_states> &>(),
	 std::declval<const system_t&>(),
	 std::declval<scalar_t>(),
	 std::declval<scalar_t>(),
	 std::declval<::pressio::ode::types::step_t>()
	 )
       )
      >::value
    >
  > : std::true_type{};
/////////////////////////

template<
  typename T, typename state_t, typename residual_t, typename system_t, typename scalar_t,
  typename enable = void
  >
struct residual_policy_callable_with_seven_args : std::false_type{};

template<
  typename T, typename state_t, typename residual_t, typename system_t, typename scalar_t
  >
struct residual_policy_callable_with_seven_args<
  T, state_t, residual_t, system_t, scalar_t,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T>().operator()
       ( std::declval<const state_t &>(),
	 std::declval<residual_t &>(),
	 std::declval<const std::array<state_t, T::num_aux_states> &>(),
	 std::declval<const system_t&>(),
	 std::declval<scalar_t>(),
	 std::declval<scalar_t>(),
	 std::declval<::pressio::ode::types::step_t>()
	 )
       )
      >::value
    >
  > : std::true_type{};
/////////////////////////

template<
  typename T,
  typename state_t,
  typename residual_t,
  typename system_t,
  typename scalar_t,
  typename enable = void
  >
struct is_legitimate_residual_policy_for_implicit_arbitrary_stepper
{
  static constexpr bool c1 = ::pressio::ode::meta::is_legitimate_implicit_state_type<state_t>::value;
  static constexpr bool c2 = ::pressio::ode::meta::is_legitimate_implicit_residual_type<residual_t>::value;

  static constexpr bool c3 = ::pressio::ode::meta::has_stepper_order_static_member<T>::value;
  static_assert( c3, "The residual policy you are trying to pass to \
arbitrary implicit stepper is missing a static member: stepper_order=... \
to set the order of the stepper.");

  static constexpr bool c4 = ::pressio::ode::meta::has_num_aux_states_static_member<T>::value;
  static_assert( c4, "The residual policy you are trying to pass to \
arbitrary implicit stepper is missing a static member: num_aux_states=... \
setting the number of auxiliary states I need.");
  static constexpr types::stepper_n_states_t numAuxStates = T::num_aux_states;

  static constexpr bool c5 = residual_policy_callable_with_six_args<
    T, state_t, residual_t, system_t, scalar_t>::value;

   static constexpr bool c6 = residual_policy_callable_with_seven_args<
    T, state_t, residual_t, system_t, scalar_t>::value;

  using value_type = bool;
  static constexpr value_type value = c1 && c2 && c3 && c4 && c5 && c6;
  using type = std::integral_constant<value_type, value>;
};

}}} // namespace pressio::ode::meta
#endif
