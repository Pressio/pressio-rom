/*
//@HEADER
// ************************************************************************
//
// ode_predicates.hpp
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

#ifndef PRESSIOROM_ODE_CONCEPTS_ODE_PREDICATES_HPP_
#define PRESSIOROM_ODE_CONCEPTS_ODE_PREDICATES_HPP_

#include <optional>

namespace pressio{ namespace ode{

template <class T, class StateType, class = void>
struct has_const_create_state_method_return_result : std::false_type{};

template <class T, class StateType>
struct has_const_create_state_method_return_result<
  T, StateType,
  std::enable_if_t<
    std::is_same<
      StateType,
      decltype(std::declval<T const>().createState())
      >::value
    >
  > : std::true_type{};


template <class T, class RhsType, class = void>
struct has_const_create_rhs_method_return_result
  : std::false_type{};

template <class T, class RhsType>
struct has_const_create_rhs_method_return_result<
  T, RhsType,
  std::enable_if_t<
    !std::is_void<RhsType>::value and
    std::is_same<
      RhsType,
      decltype(
	       std::declval<T const>().createRhs()
	       )
      >::value
    >
  > : std::true_type{};


template <class T, class MMType, class = void>
struct has_const_create_mass_matrix_method_return_result
  : std::false_type{};

template <class T, class MMType>
struct has_const_create_mass_matrix_method_return_result<
  T, MMType,
  std::enable_if_t<
    !std::is_void<MMType>::value and
    std::is_same<
      MMType,
      decltype(
	       std::declval<T const>().createMassMatrix()
	       )
      >::value
    >
  > : std::true_type{};


template <class T, class JacobianType, class = void>
struct has_const_create_jacobian_method_return_result
  : std::false_type{};

template <class T, class JacobianType>
struct has_const_create_jacobian_method_return_result<
  T, JacobianType,
  std::enable_if_t<
    !std::is_void<JacobianType>::value and
    std::is_same<
      JacobianType,
      decltype(
         std::declval<T const>().createJacobian()
         )
      >::value
    >
  > : std::true_type{};


template <class T, class ResultType, class = void>
struct has_const_create_discrete_residual_method_return_result
  : std::false_type{};

template <class T, class ResultType>
struct has_const_create_discrete_residual_method_return_result<
  T, ResultType,
  std::enable_if_t<
    !std::is_void<ResultType>::value and
    std::is_same<
      ResultType,
      decltype
      (
       std::declval<T const>().createDiscreteResidual()
       )
      >::value
    >
  > : std::true_type{};


template <class T, class JacobianType, class = void>
struct has_const_create_discrete_jacobian_method_return_result
  : std::false_type{};

template <class T, class JacobianType>
struct has_const_create_discrete_jacobian_method_return_result<
  T, JacobianType,
  std::enable_if_t<
    !std::is_void<JacobianType>::value and
    std::is_same<
      JacobianType,
      decltype(
         std::declval<T const>().createDiscreteJacobian()
         )
      >::value
    >
  > : std::true_type{};


template <
  class T, int numStates, class StepType, class IndVarType, class StateType,
  class = void
  >
struct has_const_pre_step_hook_method
  : std::false_type{};

template <class T, class StepType, class IndVarType, class StateType>
struct has_const_pre_step_hook_method<
  T, 2, StepType, IndVarType, StateType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const>().preStepHook
       (
	std::declval<StepType const &>(),
	std::declval<IndVarType const &>(),
	std::declval<IndVarType const &>(),
	std::declval<StateType const&>(),
	std::declval<StateType const&>()
	)
       )
      >::value
    >
  > : std::true_type{};



//
// has_const_discrete_residual_jacobian_method
//

template <
  class T, int n,
  class StepType, class IndVarType, class state_t,
  class ResidualType,
  class JacobianType,
  class = void
  >
struct has_const_discrete_residual_jacobian_method
  : std::false_type{};

template <
  class T, class StepType, class IndVarType, class state_t,
  class ResidualType, class JacobianType>
struct has_const_discrete_residual_jacobian_method<
  T, 1, StepType, IndVarType, state_t, ResidualType, JacobianType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const>().discreteResidualAndJacobian
       (
	std::declval<StepType const &>(),
	std::declval<IndVarType const &>(),
	std::declval<IndVarType const &>(),
	std::declval<ResidualType &>(),
  std::declval< std::optional<JacobianType*> >(),
	std::declval<state_t const&>()
	)
       )
      >::value
    >
  > : std::true_type{};


template <
  class T, class StepType, class IndVarType, class state_t,
  class ResidualType, class JacobianType>
struct has_const_discrete_residual_jacobian_method<
  T, 2, StepType, IndVarType, state_t, ResidualType, JacobianType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const>().discreteResidualAndJacobian
       (
         std::declval<StepType const &>(),
	 std::declval<IndVarType const &>(),
	 std::declval<IndVarType const &>(),
	 std::declval<ResidualType &>(),
   std::declval< std::optional<JacobianType*> >(),
	 std::declval<state_t const&>(),
	 std::declval<state_t const&>()
	)
       )
      >::value
    >
  > : std::true_type{};

template <
  class T, class StepType, class IndVarType, class state_t,
  class ResidualType, class JacobianType>
struct has_const_discrete_residual_jacobian_method<
  T, 3, StepType, IndVarType, state_t, ResidualType, JacobianType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const>().discreteResidualAndJacobian
       (
	std::declval<StepType const &>(),
	std::declval<IndVarType const &>(),
	std::declval<IndVarType const &>(),
	std::declval<ResidualType &>(),
  std::declval< std::optional<JacobianType*> >(),
	std::declval<state_t const&>(),
	std::declval<state_t const&>(),
	std::declval<state_t const&>()
	)
       )
      >::value
    >
  > : std::true_type{};


template <
  class T, class StepType, class IndVarType, class state_t,
    class ResidualType, class JacobianType>
struct has_const_discrete_residual_jacobian_method<
  T, 4, StepType, IndVarType, state_t, ResidualType, JacobianType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const>().discreteResidualAndJacobian
       (
	std::declval<StepType const &>(),
	std::declval<IndVarType const &>(),
	std::declval<IndVarType const &>(),
	std::declval<ResidualType &>(),
  std::declval< std::optional<JacobianType*> >(),
	std::declval<state_t const&>(),
	std::declval<state_t const&>(),
	std::declval<state_t const&>(),
	std::declval<state_t const&>()
	)
       )
      >::value
    >
  > : std::true_type{};


}}
#endif  // PRESSIOROM_ODE_CONCEPTS_ODE_PREDICATES_FOR_SYSTEM_HPP_
