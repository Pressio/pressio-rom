/*
//@HEADER
// ************************************************************************
//
// ode_advance_fncs.hpp
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

namespace pressio{ namespace ode{

// Convenience functions for policies
template<class Time>
inline StepsFixed<Time, policy::fixed_dt<Time>>
steps_fixed_dt(Time t0, StepCount n, Time dt){
  return { t0, n, policy::fixed_dt<Time>{dt} };
}

template<class Time, class StepSizer>
inline StepsFixed<Time, StepSizer>
steps(Time t0, StepCount n, StepSizer sizer){
  return { t0, n, std::move(sizer) };
}

template<class Time, class StepSizer>
inline ToTime<Time, StepSizer>
to_time(Time t0, Time tf, StepSizer sizer){
  static_assert(StepSizePolicy<StepSizer, Time>::value, "invalid step sizer policy");
  return { t0, tf, std::move(sizer)};
}


/* ============================================================================

advance (no solver)

High-level
----------
Advance an ODE state in time using a user-provided "Stepper" that does *not*
require a solver to perform a step.
The time marching is driven by a "Policy" object that decides:
  - when to start (start_time()),
  - when to stop (keep_going(time, step_count)),
  - the step size at each iteration (set_dt(step_start_at, step_count, dt)).

Optionally, callers may pass up to two trailing "observer" callbacks:
  - StateObserver   — called at t0 (before stepping) and after each successful step
  - RhsObserver     — optionally forwarded to the stepper iff the stepper supports it

Accepted observer forms:
  - none, (StateObserver), (StateObserver, RhsObserver)

Notes
---------------------------
1) Initialization:
     Time t = policy.start_time();
     state_obs(0, t, state);          // observe the initial condition

2) Each iteration (for step count k):
     policy.set_dt(t, k, dt);         // policy decides dt for the step starting at t
     stepper(..., t, k, dt [, rhs]);  // perform one explicit step
     t += dt;                          // advance time
     state_obs(k, t, state);          // observe the state *at the end* of step k
     ++k;

- RhsObserver is only passed to the stepper when BOTH of these hold:
    - a non-trivial rhs observer was provided, and
    - the stepper has an overload that accepts it.
- StateObserver position is intentional: it’s called once at t0 (k=0), and then
  after each step using the step count k that was used for that step. This
  makes the (k, t) pair consistent with the completed update.

Ownership & lifetime
--------------------
- Observers passed as lvalues are stored as std::reference_wrapper (caller must keep them alive).
- Observers passed as rvalues are moved and stored by value (owned here).

============================================================================ */

template<class State, class Stepper, class Policy, class... ObserversOrEmpty>
std::enable_if_t< PRESSIO_VALUE_OF(StepperWithoutSolver<Stepper>) >
advance(Stepper& stepper,
	State& state,
	const Policy& policy,
	ObserversOrEmpty&&... obs_tail)
{
  /* Deduce the time type from the policy. Validate the policy’s interface. */
  using Time = decltype(policy.start_time());

  static_assert(std::is_same_v<
		decltype(policy.keep_going(std::declval<Time>(),
					   std::declval<StepCount>())),
		bool>,
                "policy must expose keep_going(Time,StepCount)->bool");

  static_assert(std::is_void_v<
		decltype(policy.set_dt(std::declval<StepStartAt<Time>>(),
				       std::declval<StepCount>(),
				       std::declval<StepSize<Time>&>()))
		>,
                "policy must expose set_dt(Time,StepCount,StepSize<>dt)");

  static_assert(sizeof...(ObserversOrEmpty) <= 2,
    "Only nothing or observers are allowed here: none, (StateObserver), or (StateObserver,RhsObserver).");

  // Parse and hold observers with the correct ownership semantics.
  auto obs = parse_observers<Time, State>(std::forward<ObserversOrEmpty>(obs_tail)...);
  auto&& state_obs = obs.state;
  auto&& rhs_obs   = obs.rhs;
  // Detect whether a non-trivial rhs observer was provided (vs. the no-op).
  constexpr bool has_non_trivial_rhs_obs = !std::is_same< std::decay_t<decltype(rhs_obs)>, rhs_noop_t>::value;

  // Initialize time and notify the state observer about the initial condition.
  Time t = policy.start_time();
  state_obs(StepCount{0}, t, state);

  // Main advance loop. k is the step counter; its initial value is defined elsewhere
  // (e.g., first_step_value) to match the project’s global step indexing convention.
  ::pressio::ode::StepSize<Time> dt = {};
  StepCount k{first_step_value};

  while (policy.keep_going(t, k)) {
    // Form the strongly-typed "step starts at t" tag and query the step size.
    StepStartAt<Time> stepStartTime(t);
    policy.set_dt(stepStartTime, k, dt);
    impl::print_step_and_current_time(k.get(), t, dt.get());

    /* Dispatch to the stepper. If an rhs observer was provided AND the stepper
       supports an overload that takes it, call that. Otherwise call the simpler
       overload. If neither is available, trigger a compile-time error with a
       clear diagnostic. */
    if constexpr (has_non_trivial_rhs_obs &&
		  has_explicit_step_with_rhs<Stepper, State, Time, Time, decltype(rhs_obs)>::value) {
      stepper(state, stepStartTime, k, dt, rhs_obs);
    }
    else if constexpr (has_explicit_step<Stepper, State, Time, Time>::value) {
      stepper(state, stepStartTime, k, dt);
    }
    else{
      static_assert([]{return false;}(),
		    "Stepper must provide operator(state,t,dt,k) or operator(state,t,dt,k,rhsObserver).");
    }

    // Advance physical time by the step size just used
    t += dt.get();

    //Observe the state *after* completing step k. The placement here is
    //intentional to make (k, t) correspond to the just-finished step.
    state_obs(k, t, state);

    // move on to next step
    ++k;
  }
}

/* ============================================================================

advance (with solver)

High-level
----------
Advance an ODE state in time using a user-provided "Stepper" that needs a
solver to perform one step.
The time marching is driven by a "Policy" that decides:
  - when to start (start_time()),
  - when to stop (keep_going(time, step_count)),
  - the step size at each iteration (set_dt(step_start_at, step_count, dt)).

In addition, the call accepts:
  - a Solver& object, forwarded to the Stepper at every step;
  - zero or more extra solver arguments (packed and reused each step);
  - up to two trailing observer callbacks parsed from the tail:
      - StateObserver   — called at t0 and after each successful step
      - RhsObserver     — optionally forwarded to the Stepper if supported

Accepted observer forms in the tail:
  - none, (StateObserver), (StateObserver, RhsObserver)

Extra solver args
-----------------
- Any additional arguments intended for the solver/stepper are collected into a
  tuple and reused each iteration.
- They are *passed as lvalues* each time. If an rvalue is provided,
  it will be stored in the tuple by value and then used as an lvalue
  on every step; it is not "consumed" per step.

Overload selection for rhs observer
-----------------------------------
- If a non-trivial RhsObserver is present *and* the Stepper overload that accepts
  it is available, that overload is called. Otherwise, the simpler overload is used.
- Compile-time detection uses a trait on the exact call signature built from the
  current types (Stepper, State, Time, Solver, extra args..., rhs_obs).

Ownership & lifetime
--------------------
- Observers passed as lvalues are stored as std::reference_wrapper (caller keeps alive).
- Observers passed as rvalues are moved and stored by value (owned here).
- Extra solver args are stored in a tuple and reused each iteration; if they
  capture references, ensure those referents outlive the loop.

============================================================================ */

template<class State, class Stepper, class Policy, class Solver, class... Rest>
std::enable_if_t< !PRESSIO_VALUE_OF(StepperWithoutSolver<Stepper>) >
advance(Stepper& stepper,
	State& state,
	const Policy& policy,
	Solver& solver,
	Rest&&... rest)
{
  using Time = decltype(policy.start_time());

  static_assert(std::is_same_v<
		decltype(policy.keep_going(std::declval<Time>(),
					   std::declval<StepCount>())),
		bool>,
                "policy must expose keep_going(Time,StepCount)->bool");

  static_assert(std::is_void_v<
		decltype(policy.set_dt(std::declval<StepStartAt<Time>>(),
				       std::declval<StepCount>(),
				       std::declval<StepSize<Time>&>()))
		>,
                "policy must expose set_dt(Time,StepCount,StepSize<>dt)");

  auto split = split_solver_and_observers<Time, State>(solver, std::forward<Rest>(rest)...);
  auto& solver_args = split.first;
  auto  obs         = std::move(split.second);

  auto&& state_obs = obs.state;
  auto&& rhs_obs   = obs.rhs;
  constexpr bool has_non_trivial_rhs_obs = !std::is_same< std::decay_t<decltype(rhs_obs)>, rhs_noop_t>::value;

  Time t = policy.start_time();
  state_obs(StepCount{0}, t, state);

  ::pressio::ode::StepSize<Time> dt = {};
  StepCount k{first_step_value};

  while (policy.keep_going(t, k)) {
    policy.set_dt(StepStartAt<Time>(t), k, dt);
    impl::print_step_and_current_time(k.get(), t, dt.get());

    // Build a small dispatcher that calls the appropriate stepper overload.
    // The dispatcher receives the unpacked solver_args as lvalues each time.
    auto call = [&](auto&... args){
      // Can we call overload WITH rhs_obs as the last argument?
      using with_rhs = pressio::ode::detail::is_detected<
	void, // dummy to pick the specialization
	pressio::ode::detail::do_step_with_solver_expr,
	Stepper, State, Time, Time, Solver, decltype(args)..., decltype(rhs_obs)>;

      if constexpr (has_non_trivial_rhs_obs && with_rhs::value) {
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
     !PRESSIO_VALUE_OF(StepperWithoutSolver<StepperType>)
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
     !PRESSIO_VALUE_OF(StepperWithoutSolver<StepperType>)
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
