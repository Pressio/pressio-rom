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

#ifndef PRESSIOROM_ODE_IMPL_ODE_ADVANCE_N_STEPS_HPP_
#define PRESSIOROM_ODE_IMPL_ODE_ADVANCE_N_STEPS_HPP_

#include "ode_advance_printing_helpers.hpp"

namespace pressio{ namespace ode{ namespace impl{

template <
  class StepperType,
  class IndVarType,
  class dt_policy,
  class StateType,
  class ObserverType,
  class ... Args
  >
void advance_n_steps_with_dt_policy(StepperType & stepper,
				    ::pressio::ode::StepCount numSteps,
				    const IndVarType & start_val,
				    StateType & odeState,
				    dt_policy && dtManager,
				    ObserverType && observer,
				    Args && ... args)
{

  using step_t = typename ::pressio::ode::StepCount::value_type;
  IndVarType time = start_val;
  observer(StepCount{0}, start_val, odeState);

  // default construct
  ::pressio::ode::StepSize<IndVarType> dt;
  step_t step = ::pressio::ode::first_step_value;
  PRESSIOLOG_DEBUG("impl: advance_n_steps_with_dt_policy");
  for( ; step <= numSteps.get(); ++step)
    {
      const auto stepWrap = ::pressio::ode::StepCount(step);

      // call the dt manager to set the dt to use
      dtManager(stepWrap,
		::pressio::ode::StepStartAt<IndVarType>(time),
		dt);
      print_step_and_current_time(step, time, dt.get());

      stepper(odeState,
	      ::pressio::ode::StepStartAt<IndVarType>(time),
	      stepWrap, dt,
	      std::forward<Args>(args)...);

      time += dt.get();
      observer(::pressio::ode::StepCount(step), time, odeState);
    }
}


template <
  class StepperType,
  class IndVarType,
  class StateType,
  class ObserverType,
  class ... Args
  >
void advance_n_steps_with_fixed_dt(StepperType & stepper,
				   const ::pressio::ode::StepCount & numSteps,
				   const IndVarType & start_val,
				   const IndVarType & step_size,
				   StateType & odeState,
				   ObserverType && observer,
				   Args && ... args)
{


  const auto dtSetter =
    [sz = step_size](pressio::ode::StepCount /*currStep*/,
		     pressio::ode::StepStartAt<IndVarType> /*currTime*/,
		     pressio::ode::StepSize<IndVarType> & dt){
      dt = sz;
    };

  advance_n_steps_with_dt_policy(stepper, numSteps,
				 start_val, odeState, dtSetter,
				 std::forward<ObserverType>(observer),
				 std::forward<Args>(args)...);
}

}}}//end namespace pressio::ode::impl
#endif  // PRESSIOROM_ODE_IMPL_ODE_ADVANCE_N_STEPS_HPP_
