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

//
// const dt
//
template<class StepperType, class StateType, class IndVarType>
std::enable_if_t< ExplicitStepper<StepperType>::value >
advance_n_steps(StepperType & stepper,
		StateType & state,
		const IndVarType & startVal,
		const IndVarType & stepSize,
		StepCount numSteps)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, startVal);
  using observer_t = impl::NoOpStateObserver<IndVarType, StateType>;
  impl::advance_n_steps_with_fixed_dt(stepper, numSteps, startVal,
				      stepSize, state,
				      observer_t());
}

//
// dt policy provided
//
template<
  class StepperType,
  class StateType,
  class StepSizePolicyType,
  class IndVarType
  >
std::enable_if_t<
  ExplicitStepper<StepperType>::value
  && StepSizePolicy<StepSizePolicyType &&, IndVarType>::value
>
advance_n_steps(StepperType & stepper,
		StateType & state,
		const IndVarType & startVal,
		StepSizePolicyType && stepSizePolicy,
		StepCount numSteps)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, startVal);
  using observer_t = impl::NoOpStateObserver<IndVarType, StateType>;
  impl::advance_n_steps_with_dt_policy(stepper, numSteps, startVal, state,
				       std::forward<StepSizePolicyType>(stepSizePolicy),
				       observer_t());
}

//
// const dt and observer
//
template<
  class StepperType,
  class StateType,
  class ObserverType,
  class IndVarType
  >
std::enable_if_t<
  ExplicitStepper<StepperType>::value
  && StateObserver<ObserverType &&, IndVarType, StateType>::value
>
advance_n_steps(StepperType & stepper,
		StateType & state,
		const IndVarType & startVal,
		const IndVarType & stepSize,
		StepCount numSteps,
		ObserverType && observer)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, startVal);
  impl::advance_n_steps_with_fixed_dt(stepper, numSteps, startVal,
				      stepSize, state,
				      std::forward<ObserverType>(observer));
}

//
// dt policy provided and observer
//
template<
  class StepperType,
  class StateType,
  class ObserverType,
  class StepSizePolicyType,
  class IndVarType>
std::enable_if_t<
     ExplicitStepper<StepperType>::value
  && StepSizePolicy<StepSizePolicyType &&, IndVarType>::value
  && StateObserver<ObserverType &&, IndVarType, StateType>::value
>
advance_n_steps(StepperType & stepper,
		StateType & state,
		const IndVarType & startVal,
		StepSizePolicyType && stepSizePolicy,
		StepCount numSteps,
		ObserverType && observer)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, startVal);
  impl::advance_n_steps_with_dt_policy(stepper, numSteps, startVal, state,
				       std::forward<StepSizePolicyType>(stepSizePolicy),
				       std::forward<ObserverType>(observer));
}

}} //end namespace pressio::ode
#endif  // PRESSIOROM_ODE_ODE_ADVANCE_N_STEPS_HPP_
