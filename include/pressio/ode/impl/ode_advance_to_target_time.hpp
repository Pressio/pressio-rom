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

#ifndef PRESSIOROM_ODE_IMPL_ODE_ADVANCE_TO_TARGET_TIME_HPP_
#define PRESSIOROM_ODE_IMPL_ODE_ADVANCE_TO_TARGET_TIME_HPP_

#include "ode_advance_printing_helpers.hpp"

namespace pressio{ namespace ode{ namespace impl{

template <
  bool useExtraArgs,
  class StepSizePolicyType, class IndVarType, class ...Args
  >
std::enable_if_t< useExtraArgs==true >
call_dt_policy(StepSizePolicyType && dtPolicy,
	       const StepCount & step,
	       const ::pressio::ode::StepStartAt<IndVarType> & time,
	       ::pressio::ode::StepSize<IndVarType> & dt,
	       Args && ...args)
{
  dtPolicy(step, time, dt, std::forward<Args>(args)...);
}

template<
  bool useExtraArgs,
  class StepSizePolicyType, class IndVarType, class ...Args
  >
std::enable_if_t< useExtraArgs==false >
call_dt_policy(StepSizePolicyType && dtPolicy,
	       const StepCount & step,
	       const ::pressio::ode::StepStartAt<IndVarType> time,
	       ::pressio::ode::StepSize<IndVarType> & dt,
	       Args && ...args)
{
  dtPolicy(step, time, dt);
}

template <
  bool enableTimeStepRecovery,
  class StepperType,
  class IndVarType,
  class StateType,
  class ObserverType,
  class StepSizePolicyType,
  class ... Args>
void to_target_time_with_step_size_policy(StepperType & stepper,
						const IndVarType & start_time,
						const IndVarType & final_time,
						StateType & odeState,
						StepSizePolicyType	&& dtPolicy,
						ObserverType && observer,
						Args && ... args)
{

  if (final_time < start_time){
    throw std::runtime_error("You cannot call the advancer with final time < start time.");
  }

  if (final_time == start_time){
    return;
  }

  using step_t = typename StepCount::value_type;

  IndVarType time  = start_time;

  ::pressio::ode::StepSize<IndVarType> dt{0};
  ::pressio::ode::StepSizeMinAllowedValue<IndVarType> minDt{0};
  ::pressio::ode::StepSizeScalingFactor<IndVarType> dtScalingFactor{1};

  // observe initial condition
  observer(StepCount{0}, time, odeState);

  step_t step = ::pressio::ode::first_step_value;
  PRESSIOLOG_DEBUG("impl: advance_to_target_time_with_dt_policy");
  constexpr auto eps = std::numeric_limits<IndVarType>::epsilon();
  bool condition = true;
  while (condition)
    {
      const auto stepWrap = ::pressio::ode::StepCount(step);

      PRESSIOLOG_DEBUG("callback dt policy");
      impl::call_dt_policy<enableTimeStepRecovery>(dtPolicy, stepWrap,
					   ::pressio::ode::StepStartAt<IndVarType>(time),
					   dt, minDt, dtScalingFactor);

      if (dt.get() < minDt.get()){
	throw std::runtime_error("The time step size cannot be smaller than the minimum value.");
      }

      if (enableTimeStepRecovery){
	if (dtScalingFactor.get() <= static_cast<IndVarType>(1)){
	  // need to change this to use some notion of identity
	  throw std::runtime_error("The time step size reduction factor must be > 1.");
	}
      }

      print_step_and_current_time(step, time, dt.get());

      if (enableTimeStepRecovery)
      {
	bool needStop = false;
	while(!needStop){
	  try
	  {
	    stepper(odeState,
		    ::pressio::ode::StepStartAt<IndVarType>(time),
		    stepWrap, dt,
		    std::forward<Args>(args)...);
	    needStop=true;
	  }
	  catch (::pressio::eh::TimeStepFailure const & e)
	  {
	    dt = dt.get()/dtScalingFactor.get();
	    if (dt.get() < minDt.get()){
	      throw std::runtime_error("Violation of minimum time step while trying to recover time step");
	    }

	    PRESSIOLOG_WARNING("time step={} failed, retrying with dt={}", step, dt.get());
	  }
	}
      }
      else
      {
	stepper(odeState,
		::pressio::ode::StepStartAt<IndVarType>(time),
		stepWrap, dt,
		std::forward<Args>(args)...);
      }

      time += dt.get();
      observer(::pressio::ode::StepCount(step), time, odeState);

      // use numeric limits to avoid tricky roundoff accumulation
      if ( std::abs(time - final_time) <= eps ) condition = false;

      // if we are over the final time, stop too
      if ( time > final_time ) condition = false;

      step++;
    }
}

}}}//end namespace pressio::ode::impl
#endif  // PRESSIOROM_ODE_IMPL_ODE_ADVANCE_TO_TARGET_TIME_HPP_
