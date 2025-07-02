/*
//@HEADER
// ************************************************************************
//
// ode_advance_to_target_time.hpp
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

#ifndef PRESSIOROM_ODE_ODE_ADVANCE_TO_TARGET_POINT_WITH_STEP_RECOVERY_VARIADIC_HPP_
#define PRESSIOROM_ODE_ODE_ADVANCE_TO_TARGET_POINT_WITH_STEP_RECOVERY_VARIADIC_HPP_

#include "./impl/ode_advance_noop_observer.hpp"
#include "./impl/ode_advance_to_target_time.hpp"
#include "./impl/ode_advance_mandates.hpp"

namespace pressio{ namespace ode{

template<
  class StepperType,
  class StateType,
  class StepSizePolicyType,
  class IndVarType,
  class SolverType,
  class ...SolverArgs
  >
std::enable_if_t<
     ImplicitStepper<void, StepperType, SolverType&&, SolverArgs&&...>::value
  && StepSizePolicyWithReductionScheme<StepSizePolicyType&&, IndVarType>::value
  && !StateObserver<SolverType&&, IndVarType, StateType>::value
  >
advance_to_target_point_with_step_recovery(StepperType & stepper,
					  StateType & state,
					  const IndVarType & startVal,
					  const IndVarType & finalVal,
					  StepSizePolicyType && stepSizePolicy,
					  SolverType && solver,
					  SolverArgs && ... solverArgs)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, startVal);
  using observer_t = impl::NoOpStateObserver<IndVarType, StateType>;
  impl::to_target_time_with_step_size_policy
    <true>(stepper, startVal,
	   finalVal, state,
	   std::forward<StepSizePolicyType>(stepSizePolicy),
	   observer_t(),
	   std::forward<SolverType>(solver),
	   std::forward<SolverArgs>(solverArgs)...);
}

template<
  class StepperType,
  class StateType,
  class StepSizePolicyType,
  class ObserverType,
  class IndVarType,
  class SolverType,
  class ...SolverArgs
  >
std::enable_if_t<
     ImplicitStepper<void, StepperType, SolverType&&, SolverArgs&&...>::value
  && StepSizePolicyWithReductionScheme<StepSizePolicyType&&, IndVarType>::value
  && StateObserver<ObserverType&&, IndVarType, StateType>::value
  >
advance_to_target_point_with_step_recovery(StepperType & stepper,
					  StateType & state,
					  const IndVarType & startVal,
					  const IndVarType & finalVal,
					  StepSizePolicyType && stepSizePolicy,
					  ObserverType && observer,
					  SolverType && solver,
					  SolverArgs && ... solverArgs)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, startVal);
  impl::to_target_time_with_step_size_policy
    <true>(stepper, startVal,
	   finalVal, state,
	   std::forward<StepSizePolicyType>(stepSizePolicy),
	   std::forward<ObserverType>(observer),
	   std::forward<SolverType>(solver),
	   std::forward<SolverArgs>(solverArgs)...);
}

}}//end namespace pressio::ode
#endif  // PRESSIOROM_ODE_ODE_ADVANCE_TO_TARGET_POINT_WITH_STEP_RECOVERY_VARIADIC_HPP_
