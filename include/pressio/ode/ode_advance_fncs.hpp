/*
//@HEADER
// ************************************************************************
//
// ode_advance_variadic.hpp
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

#ifndef PRESSIOROM_ODE_ODE_ADVANCE_FNCS_HPP_
#define PRESSIOROM_ODE_ODE_ADVANCE_FNCS_HPP_

#include "./impl/ode_advance_impl.hpp"

#include <tuple>
#include <type_traits>
#include <utility>
#include <functional>

namespace pressio{ namespace ode{

// Convenience functions for policies
template<class Time>
inline StepsFixed<Time, policy::fixed_dt<Time>>
steps_fixed_dt(Time t0, StepCount n, Time dt){ return { t0, n, policy::fixed_dt<Time>{dt} }; }

template<class Time, class StepSizer>
inline StepsFixed<Time, StepSizer>
steps(Time t0, StepCount n, StepSizer sizer){ return { t0, n, std::move(sizer) }; }

template<class Time, class StepSizer>
inline ToTime<Time, StepSizer>
to_time(Time t0, Time tf, StepSizer sizer){ return { t0, tf, std::move(sizer)}; }


//======================== advance (no solver) ====================
template<class State, class Stepper, class Policy, class... Obs>
std::enable_if_t< StepperWithoutSolver<Stepper>::value >
advance(Stepper& stepper,
	State& state,
	const Policy& bp,
	Obs&&... obs_tail)
{
  using Time = decltype(bp.start_time());
  static_assert(std::is_same_v<decltype(bp.keep_going(std::declval<Time>(), std::declval<StepCount>())), bool>,
                "boundsPolicy must expose keep_going(Time,StepCount)->bool");
  static_assert(std::is_void_v<decltype(bp.next_dt(std::declval<StepStartAt<Time>>(), std::declval<StepCount>(), std::declval<StepSize<Time>&>()))>,
                "boundsPolicy must expose next_dt(t,k,dt)");

  static_assert(sizeof...(Obs) <= 2,
    "Only observers are allowed here: none, (StateObserver), or (StateObserver,RhsObserver).");

  auto obs = parse_observers<Time, State>(std::forward<Obs>(obs_tail)...);
  auto&& state_obs = obs.state;
  auto&& rhs_obs   = obs.rhs;
  constexpr bool have_rhs = !std::is_same< std::decay_t<decltype(rhs_obs)>, rhs_noop_t>::value;

  Time t = bp.start_time();
  state_obs(StepCount{0}, t, state);

  ::pressio::ode::StepSize<Time> dt;
  StepCount k{first_step_value};
  while (bp.keep_going(t, k)) {
    auto stepStartVal = StepStartAt<Time>(t);
    bp.next_dt(stepStartVal, k, dt);
    impl::print_step_and_current_time(k.get(), t, dt.get());

    // Prefer passing rhs_obs into the stepper if supported and provided
    if constexpr (have_rhs &&
		  has_explicit_step_with_rhs<Stepper, State, Time, Time, decltype(rhs_obs)>::value) {
      stepper(state, stepStartVal, k, dt, rhs_obs);
    }
    else if constexpr (has_explicit_step<Stepper, State, Time, Time>::value) {
      stepper(state, stepStartVal, k, dt);
    }
    else{
      static_assert([]{return false;}(),
		    "Stepper must provide operator(state,t,dt,k) or operator(state,t,dt,k,rhsObserver).");
    }

    t += dt.get();
    state_obs(k, t, state);
    ++k;
  }
}

//====================== advance (with solver) ====================
template<class State, class Stepper, class Policy, class Solver, class... Rest>
std::enable_if_t< !StepperWithoutSolver<Stepper>::value >
advance(Stepper& stepper,
	State& state,
	const Policy& bp,
	Solver& solver,
	Rest&&... rest)
{
  using Time = decltype(bp.start_time());
  auto split = split_solver_and_observers<Time, State>(solver, std::forward<Rest>(rest)...);
  auto& solver_args = split.first;
  auto  obs         = std::move(split.second);

  auto&& state_obs = obs.state;
  auto&& rhs_obs   = obs.rhs;
  constexpr bool have_rhs = !std::is_same< std::decay_t<decltype(rhs_obs)>, rhs_noop_t>::value;

  Time t = bp.start_time();
  state_obs(StepCount{0}, t, state);

  ::pressio::ode::StepSize<Time> dt;
  StepCount k{first_step_value};
  while (bp.keep_going(t, k)) {
    bp.next_dt(StepStartAt<Time>(t), k, dt);
    impl::print_step_and_current_time(k.get(), t, dt.get());

    auto call = [&](auto&... args){
      // Can we call overload WITH rhs_obs as the last argument?
      using with_rhs = pressio::ode::detail::is_detected<
	void, // dummy to pick the specialization
	pressio::ode::detail::do_step_with_solver_expr,
	Stepper, State, Time, Time, Solver, decltype(args)..., decltype(rhs_obs)>;

      if constexpr (have_rhs && with_rhs::value) {
	stepper(state, StepStartAt(t), k, dt, solver, args..., rhs_obs);
      } else{
	stepper(state, StepStartAt(t), k, dt, solver, args...);
      }
    };
    std::apply(call, solver_args);

    t += dt.get();
    state_obs(k, t, state);
    ++k;
  }
}


template<
  class StepperType,
  class StateType,
  class StepSizePolicyType,
  class IndVarType,
  class SolverType,
  class ...SolverArgs
  >
std::enable_if_t<
     !StepperWithoutSolver<StepperType>::value
  && StepSizePolicyWithReductionScheme<StepSizePolicyType&&, IndVarType>::value
  && !StateObserver<SolverType&&, IndVarType, StateType>::value
  >
advance_with_step_recovery(StepperType & stepper,
			   StateType & state,
			   const IndVarType & startVal,
			   const IndVarType & finalVal,
			   StepSizePolicyType && stepSizePolicy,
			   SolverType && solver,
			   SolverArgs && ... solverArgs)
{

  impl::mandate_on_ind_var_and_state_types(stepper, state, startVal);
  using observer_t = state_noop_t;
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
     !StepperWithoutSolver<StepperType>::value
  && StepSizePolicyWithReductionScheme<StepSizePolicyType&&, IndVarType>::value
  && StateObserver<ObserverType&&, IndVarType, StateType>::value
  >
advance_with_step_recovery(StepperType & stepper,
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

}} //end namespace pressio::ode
#endif  // PRESSIOROM_ODE_ODE_ADVANCE_N_STEPS_HPP_
