/*
//@HEADER
// ************************************************************************
//
// ode_has_const_discrete_residual_jacobian_action_method.hpp
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

#ifndef PRESSIOROM_ROM_IMPL_ODE_HAS_CONST_DISCRETE_RESIDUAL_JACOBIAN_ACTION_METHOD_HPP_
#define PRESSIOROM_ROM_IMPL_ODE_HAS_CONST_DISCRETE_RESIDUAL_JACOBIAN_ACTION_METHOD_HPP_

namespace pressio{ namespace rom{

template <
  class T,
  int n,
  class StepType,
  class IndVarType,
  class StateType,
  class ResidualType,
  class ManifoldJacobian,
  class JacobianType,
  class = void
  >
struct has_const_discrete_residual_jacobian_action_method : std::false_type{};

template <
  class T, class StepType, class IndVarType, class StateType,
  class ResidualType, class ManifoldJacobian, class JacobianType>
struct has_const_discrete_residual_jacobian_action_method<
  T, 1, StepType, IndVarType, StateType, ResidualType, ManifoldJacobian, JacobianType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const>().discreteTimeResidualAndJacobianAction
       (
	  std::declval<StepType const &>(),
	  std::declval<IndVarType const &>(),
	  std::declval<IndVarType const &>(),
	  std::declval<ResidualType &>(),
	  std::declval<ManifoldJacobian const &>(),
    std::declval< std::optional<JacobianType*> >(),
	  std::declval<StateType const&>()
	)
       )
      >::value
    >
  > : std::true_type{};

template <
  class T, class StepType, class IndVarType, class StateType,
  class ResidualType, class ManifoldJacobian, class JacobianType>
struct has_const_discrete_residual_jacobian_action_method<
  T, 2, StepType, IndVarType, StateType, ResidualType, ManifoldJacobian, JacobianType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
      std::declval<T const>().discreteTimeResidualAndJacobianAction
      (
       std::declval<StepType const &>(),
       std::declval<IndVarType const &>(),
       std::declval<IndVarType const &>(),
       std::declval<ResidualType &>(),
       std::declval<ManifoldJacobian const &>(),
       std::declval< std::optional<JacobianType*> >(),
       std::declval<StateType const&>(),
       std::declval<StateType const&>()
       )
       )
      >::value
  >
  > : std::true_type{};

template <
  class T, class StepType, class IndVarType, class StateType,
  class ResidualType, class ManifoldJacobian, class JacobianType>
struct has_const_discrete_residual_jacobian_action_method<
  T, 3, StepType, IndVarType, StateType, ResidualType, ManifoldJacobian, JacobianType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
      std::declval<T const>().discreteTimeResidualAndJacobianAction
      (
       std::declval<StepType const &>(),
       std::declval<IndVarType const &>(),
       std::declval<IndVarType const &>(),
       std::declval<ResidualType &>(),
       std::declval<ManifoldJacobian const &>(),
       std::declval< std::optional<JacobianType*> >(),
       std::declval<StateType const&>(),
       std::declval<StateType const&>(),
       std::declval<StateType const&>()
       )
       )
  >::value
  >
  > : std::true_type{};

}}
#endif  // PRESSIOROM_ROM_IMPL_ODE_HAS_CONST_DISCRETE_RESIDUAL_JACOBIAN_ACTION_METHOD_HPP_
