/*
//@HEADER
// ************************************************************************
//
// ode_others.hpp
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

#ifndef PRESSIOROM_ODE_CONCEPTS_ODE_OTHERS_HPP_
#define PRESSIOROM_ODE_CONCEPTS_ODE_OTHERS_HPP_

namespace pressio{ namespace ode{

template <class T, class = void>
struct ExplicitStepper : std::false_type{};

template <class T>
struct ExplicitStepper<
  T,
  std::enable_if_t<
    ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval< typename T::state_type & >(),
	std::declval< ::pressio::ode::StepStartAt<typename T::independent_variable_type> >(),
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepSize<typename T::independent_variable_type> >()
	)
       )
      >::value
    //&& impl::stepper_accepting_lvalue_state<T>::value
    >
  > : std::true_type{};

template <class T, class AuxT, class ...Args>
struct ImplicitStepper : std::false_type{};

template <class T, class AuxT, class ...Args>
struct ImplicitStepper<
  std::enable_if_t<
    ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval<typename T::state_type & >(),
	std::declval<::pressio::ode::StepStartAt<typename T::independent_variable_type> >(),
	std::declval<::pressio::ode::StepCount >(),
	std::declval<::pressio::ode::StepSize<typename T::independent_variable_type> >(),
	std::declval< AuxT >(), std::declval<Args>()...
	)
       )
      >::value
    //&& impl::variadic_stepper_accepting_lvalue_state<void, T, AuxT, Args...>::value
    >,
  T, AuxT, Args...
  > : std::true_type{};


template <class T, class IndVarType, class StateType, class enable = void>
struct StateObserver : std::false_type{};

template <class T, class IndVarType, class StateType>
struct StateObserver<
  T, IndVarType, StateType,
  std::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T>().operator()
	       (
		std::declval< ::pressio::ode::StepCount >(),
		std::declval< IndVarType >(),
		std::declval<StateType const &>()
		)
	       )
      >::value
    >
  > : std::true_type{};

template <class T, class IndVarType, class Enable = void>
struct StepSizePolicy : std::false_type{};

template <class T, class IndVarType>
struct StepSizePolicy<
  T, IndVarType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepStartAt<IndVarType> >(),
	std::declval< ::pressio::ode::StepSize<IndVarType> & >()
	)
       )
      >::value
    //&& impl::step_size_policy_taking_dt_by_ref<T, IndVarType>::value
    >
  > : std::true_type{};


template <class T, class IndVarType, class Enable = void>
struct StepSizePolicyWithReductionScheme : std::false_type{};

template <class T, class IndVarType>
struct StepSizePolicyWithReductionScheme<
  T, IndVarType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepStartAt<IndVarType> >(),
	std::declval< ::pressio::ode::StepSize<IndVarType> & >(),
	std::declval< ::pressio::ode::StepSizeMinAllowedValue<IndVarType> & >(),
	std::declval< ::pressio::ode::StepSizeScalingFactor<IndVarType> & >()
	)
       )
      >::value
    //&& impl::step_size_policy_with_reduc_taking_dt_by_ref<T, IndVarType>::value
    >
  > : std::true_type{};

}} // end namespace pressio::ode
#endif  // PRESSIOROM_ODE_CONCEPTS_ODE_OTHERS_HPP_
