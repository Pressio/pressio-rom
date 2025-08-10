/*
//@HEADER
// ************************************************************************
//
// ode_advance_n_steps.hpp
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

#ifndef PRESSIOROM_ODE_ODE_ADVANCE_N_STEPS_HPP_
#define PRESSIOROM_ODE_ODE_ADVANCE_N_STEPS_HPP_

#include "./impl/ode_advance_noop_observer.hpp"
#include "./impl/ode_advance_n_steps.hpp"
#include "./impl/ode_advance_mandates.hpp"

namespace pressio{ namespace ode{

template<
  class StepperType,
  class StateType,
  class IndVarType,
  class FirstArg,        // either const dt (IndVarType) or a StepSizePolicy
  class... Rest          // 0, 1, or 2: [StateObserver], [StateObserver, RhsObserver]
>
std::enable_if_t< ExplicitStepper<StepperType>::value >
advance_n_steps(StepperType & stepper,
                StateType & state,
                const IndVarType & startVal,
                FirstArg && first,
                StepCount numSteps,
                Rest&&... rest)
{
  static_assert(sizeof...(Rest) <= 2,
    "advance_n_steps: only (no observers), (StateObserver), or (StateObserver, RhsObserver) are allowed.");

  impl::mandate_on_ind_var_and_state_types(stepper, state, startVal);

  using state_noop_t = impl::NoOpStateObserver<IndVarType, StateType>;
  using rhs_noop_t   = impl::NoOpRhsObserver<IndVarType, StateType>;

  // Parse optional observers
  auto make_observers = [] (auto&&... tail)
  {
    // if the pack is empty, it means we don't want any observer
    if constexpr (sizeof...(tail) == 0) {
      return std::pair<state_noop_t, rhs_noop_t>{};
    }
    // if there is only one, then it must be a state observer
    else if constexpr (sizeof...(tail) == 1) {
      auto && s = (std::forward<decltype(tail)>(tail), ...);
      static_assert(StateObserver<mpl::remove_cvref_t<decltype(s)>, IndVarType, StateType>::value,
        "advance_n_steps: single trailing argument must be a valid StateObserver.");

      return std::pair<decltype(s), rhs_noop_t>{ std::forward<decltype(s)>(s), rhs_noop_t{} };
    }
    else { // exactly 2
      auto tup = std::forward_as_tuple(std::forward<decltype(tail)>(tail)...);
      auto && s = std::get<0>(tup);
      auto && r = std::get<1>(tup);
      static_assert(StateObserver<mpl::remove_cvref_t<decltype(s)>, IndVarType, StateType>::value,
        "advance_n_steps: first trailing argument must be a valid StateObserver.");

      static_assert(RhsObserver<mpl::remove_cvref_t<decltype(r)>, IndVarType, StateType>::value,
        "advance_n_steps: second trailing argument must be a valid RhsObserver.");
      return std::pair<decltype(s), decltype(r)>{ std::forward<decltype(s)>(s),
                                                  std::forward<decltype(r)>(r) };
    }
  };
  auto [stateObs, rhsObs] = make_observers(std::forward<Rest>(rest)...);

  // Dispatch on policy vs constant dt
  if constexpr (StepSizePolicy<FirstArg&&, IndVarType>::value) {
    impl::advance_n_steps_with_dt_policy(stepper,
                                         numSteps,
                                         startVal,
                                         state,
                                         std::forward<FirstArg>(first),
                                         std::forward<decltype(stateObs)>(stateObs),
                                         std::forward<decltype(rhsObs)>(rhsObs));
  } else {
    static_assert(std::is_convertible_v<std::decay_t<FirstArg>, IndVarType>,
      "advance_n_steps: when not passing a StepSizePolicy, `first` must be convertible to IndVarType (step size).");
    const IndVarType & stepSize = first;
    impl::advance_n_steps_with_dt_policy(stepper,
					 numSteps,
					 startVal,
					 state,
					 [sz = stepSize](StepCount /*currStep*/,
							 StepStartAt<IndVarType> /*currTime*/,
							 StepSize<IndVarType> & dt){
					   dt = sz;
					 },
					 std::forward<decltype(stateObs)>(stateObs),
					 std::forward<decltype(rhsObs)>(rhsObs));
  }
}


// //
// // const dt
// //
// template<class StepperType, class StateType, class IndVarType>
// std::enable_if_t< ExplicitStepper<StepperType>::value >
// advance_n_steps(StepperType & stepper,
// 		StateType & state,
// 		const IndVarType & startVal,
// 		const IndVarType & stepSize,
// 		StepCount numSteps)
// {

//   impl::mandate_on_ind_var_and_state_types(stepper, state, startVal);
//   using observer_t = impl::NoOpStateObserver<IndVarType, StateType>;
//   impl::advance_n_steps_with_fixed_dt(stepper, numSteps, startVal,
// 				      stepSize, state,
// 				      observer_t());
// }

// //
// // dt policy provided
// //
// template<
//   class StepperType,
//   class StateType,
//   class StepSizePolicyType,
//   class IndVarType
//   >
// std::enable_if_t<
//   ExplicitStepper<StepperType>::value
//   && StepSizePolicy<StepSizePolicyType &&, IndVarType>::value
// >
// advance_n_steps(StepperType & stepper,
// 		StateType & state,
// 		const IndVarType & startVal,
// 		StepSizePolicyType && stepSizePolicy,
// 		StepCount numSteps)
// {

//   impl::mandate_on_ind_var_and_state_types(stepper, state, startVal);
//   using observer_t = impl::NoOpStateObserver<IndVarType, StateType>;
//   impl::advance_n_steps_with_dt_policy(stepper, numSteps, startVal, state,
// 				       std::forward<StepSizePolicyType>(stepSizePolicy),
// 				       observer_t());
// }

// //
// // const dt and observer
// //
// template<
//   class StepperType,
//   class StateType,
//   class ObserverType,
//   class IndVarType
//   >
// std::enable_if_t<
//   ExplicitStepper<StepperType>::value
//   && StateObserver<ObserverType &&, IndVarType, StateType>::value
// >
// advance_n_steps(StepperType & stepper,
// 		StateType & state,
// 		const IndVarType & startVal,
// 		const IndVarType & stepSize,
// 		StepCount numSteps,
// 		ObserverType && observer)
// {

//   impl::mandate_on_ind_var_and_state_types(stepper, state, startVal);
//   impl::advance_n_steps_with_fixed_dt(stepper, numSteps, startVal,
// 				      stepSize, state,
// 				      std::forward<ObserverType>(observer));
// }

// //
// // dt policy provided and observer
// //
// template<
//   class StepperType,
//   class StateType,
//   class ObserverType,
//   class StepSizePolicyType,
//   class IndVarType>
// std::enable_if_t<
//      ExplicitStepper<StepperType>::value
//   && StepSizePolicy<StepSizePolicyType &&, IndVarType>::value
//   && StateObserver<ObserverType &&, IndVarType, StateType>::value
// >
// advance_n_steps(StepperType & stepper,
// 		StateType & state,
// 		const IndVarType & startVal,
// 		StepSizePolicyType && stepSizePolicy,
// 		StepCount numSteps,
// 		ObserverType && observer)
// {

//   impl::mandate_on_ind_var_and_state_types(stepper, state, startVal);
//   impl::advance_n_steps_with_dt_policy(stepper, numSteps, startVal, state,
// 				       std::forward<StepSizePolicyType>(stepSizePolicy),
// 				       std::forward<ObserverType>(observer));
// }

}} //end namespace pressio::ode
#endif  // PRESSIOROM_ODE_ODE_ADVANCE_N_STEPS_HPP_
