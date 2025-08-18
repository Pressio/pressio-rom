/*
//@HEADER
// ************************************************************************
//
// ode_advance_impl.hpp
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

#ifndef PRESSIOROM_ODE_IMPL_ODE_ADVANCE_IMPL_HPP_
#define PRESSIOROM_ODE_IMPL_ODE_ADVANCE_IMPL_HPP_

#include <tuple>
#include <type_traits>
#include <utility>
#include <functional>

namespace pressio{ namespace ode{ namespace impl{

template<class StepperType, class StateType, class IndVarType>
constexpr void mandate_on_ind_var_and_state_types(const StepperType & /*unused*/,
              const StateType & /*unused*/,
              const IndVarType & /*unused*/)
{
  static_assert(std::is_same<IndVarType,
    typename StepperType::independent_variable_type>::value,
    "IndVarType must be the same as StepperType::independent_variable_type");
  static_assert(std::is_same<StateType,
    typename StepperType::state_type>::value,
    "StateType must be the same as StepperType::state_type");
}

template <typename IndVarType>
void print_step_and_current_time(const typename StepCount::value_type & step,
				 const IndVarType & time,
				 const IndVarType & dt)
{
  PRESSIOLOG_DEBUG("starting timestep={} from time={} with dt={}",
		   step, time, dt);
}


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
} //end namespace impl

/* -----------------------------------------------------------------------------
   Observer parsing utilities
 ----------------------------------------------------------------------------- */

// These are default stand-ins used when the user does not supply observers.
// They accept any call signature and do nothing.
struct state_noop_t { template<class... A> void operator()(A&&...) const noexcept {} };
struct rhs_noop_t   { template<class... A> void operator()(A&&...) const noexcept {} };

/*
Storage policy for observers (obs_holder_t)

Given a deduced observer type T (which might be T&, const T&, or T&&),
choose how to *store* it inside our parsed bundle:

  • If T is an lvalue reference (T& or const T&), store std::reference_wrapper<T>.
    - This avoids copying the observer and preserves the caller's lifetime.
    - NOTE: reference_wrapper models a reference; it does *not* extend lifetime.
  • Otherwise (rvalue / prvalue), store the observer by value (decayed).
    - This takes ownership (via move) and is safe for temporaries.
*/
template<class T>
using obs_holder_t =
  std::conditional_t<
  std::is_lvalue_reference_v<T>,
  std::reference_wrapper<std::remove_reference_t<T>>,
  std::remove_reference_t<T>
  >;

/*
Helper that *constructs* an obs_holder_t from a forwarding reference.
  • If called with an lvalue, returns std::ref(v) to build a reference_wrapper.
  • If called with an rvalue, returns std::move(v) so we can store by value.

This ensures the "holder" we create matches obs_holder_t's storage policy.
*/
template<class T>
constexpr auto hold(T&& v) {
  if constexpr (std::is_lvalue_reference_v<T>) return std::ref(v);
  else                                         return std::move(v);
}

// Convenience traits that strip cv/ref before checking observer concepts
template<class Time, class State, class F>
inline constexpr bool is_state_obs_v =
  StateObserver<mpl::remove_cvref_t<F>, Time, State>::value;

template<class Time, class Rhs, class F>
inline constexpr bool is_rhs_obs_v =
  RhsObserver<mpl::remove_cvref_t<F>, Time, Rhs>::value;

// A simple bundle to return from parse_observers(...)
// parsed_obs_t holds the two observers (state and rhs) in their "holder" form.
// The types SO and RO are typically either std::reference_wrapper<...> or
// a decayed-by-value type, depending on how the user passed the observers.
template<class SO, class RO>
struct parsed_obs_t { SO state; RO rhs; };

/*
  parse_observers(...) overload set
  These constexpr overloads accept 0/1/2 trailing arguments and:

  - Validate their types at compile time via static_assert.
  - Return a parsed_obs_t with the correctly held observers.
  -For missing observers, we fill in the corresponding no-op.

IMPORTANT: The "must be placed last" rule is a *convention* enforced by how
           you integrate these helpers into your API. This overload set only
           checks count and types; your higher-level function should ensure
           it forwards the *trailing* args into parse_observers(...).
*/

// No observers provided -> use both no-ops.
template<class Time, class State>
constexpr auto parse_observers(){
  return parsed_obs_t<state_noop_t, rhs_noop_t>{ state_noop_t{}, rhs_noop_t{} };
}

// Single observer -> must be a StateObserver.
template<class Time, class State, class A1>
constexpr auto parse_observers(A1&& a1){
  static_assert(is_state_obs_v<Time, State, A1>,
    "Single observer must be a StateObserver, placed last.");
  using SOH = obs_holder_t<A1&&>;
  return parsed_obs_t<SOH, rhs_noop_t>{ hold(std::forward<A1>(a1)), rhs_noop_t{} };
}

// Two observers -> must be (StateObserver, RhsObserver) in this order.
template<class Time, class State, class A1, class A2>
constexpr auto parse_observers(A1&& a1, A2&& a2){
  static_assert(is_state_obs_v<Time, State, A1> && is_rhs_obs_v<Time, State, A2>,
    "Two observers must be (StateObserver, RhsObserver), placed last.");
  using SOH = obs_holder_t<A1&&>;
  using ROH = obs_holder_t<A2&&>;
  return parsed_obs_t<SOH, ROH>{ hold(std::forward<A1>(a1)), hold(std::forward<A2>(a2)) };
}

// ==================
//
// split solver args vs observers
//
//===============
namespace detail {

// types-only check: any observer among first K elements
template<class Time, class State, class Types, std::size_t... I>
constexpr bool any_observer_in_prefix_types(std::index_sequence<I...>) {
  return ((is_state_obs_v<Time, State, typename std::tuple_element<I, Types>::type> ||
           is_rhs_obs_v  <Time, State, typename std::tuple_element<I, Types>::type>) || ...);
}

template<std::size_t K, class Tup, std::size_t... I>
constexpr auto tuple_take_impl(Tup&& tup, std::index_sequence<I...>) {
  return std::forward_as_tuple(std::get<I>(std::forward<Tup>(tup))...);
}
template<std::size_t K, class Tup>
constexpr auto tuple_take(Tup&& tup) {
  if constexpr (K == 0) return std::tuple<>{};
  else                  return tuple_take_impl<K>(std::forward<Tup>(tup),
                                                  std::make_index_sequence<K>{});
}

} // namespace detail

template<class Time, class State, class Solver, class... Rest>
constexpr auto split_solver_and_observers(Solver&, Rest&&... rest)
{
  constexpr std::size_t N = sizeof...(Rest);
  using Types = std::tuple<std::remove_reference_t<Rest>...>;

  constexpr bool last_is_state = []() {
    if constexpr (N >= 1) {
      using L = typename std::tuple_element<N-1, Types>::type;
      return is_state_obs_v<Time, State, L>;
    } else {
      return false;
    }
  }();

  constexpr bool last_is_rhs = []() {
    if constexpr (N >= 1) {
      using L = typename std::tuple_element<N-1, Types>::type;
      return is_rhs_obs_v<Time, State, L>;
    } else {
      return false;
    }
  }();

  constexpr bool last2_is_state = []() {
    if constexpr (N >= 2) {
      using L = typename std::tuple_element<N-2, Types>::type;
      return is_state_obs_v<Time, State, L>;
    } else {
      return false;
    }
  }();

  static_assert(!(last_is_rhs && !(N >= 2 && last2_is_state)),
    "RhsObserver cannot appear without a preceding StateObserver.");

  constexpr std::size_t obs_count =
      (N >= 2 && last2_is_state && last_is_rhs) ? 2u
    : (N >= 1 && last_is_state && !last_is_rhs) ? 1u
    : 0u;

  // ensure observers are trailing only
  if constexpr (obs_count <= N) {
    constexpr std::size_t K = N - obs_count;
    static_assert(!pressio::ode::detail::any_observer_in_prefix_types<Time, State, Types>(
                    std::make_index_sequence<K>{}),
      "Observers must be the last arguments: none, (StateObserver), or (StateObserver,RhsObserver).");
  }

  // tuple of references to original args (no owning locals captured)
  auto refs = std::forward_as_tuple(std::forward<Rest>(rest)...);

  if constexpr (obs_count == 2) {
    auto solver_args = pressio::ode::detail::tuple_take<N-2>(refs);
    auto observers   = parse_observers<Time, State>(
                         std::get<N-2>(refs), std::get<N-1>(refs));
    return std::make_pair(std::move(solver_args), std::move(observers));
  } else if constexpr (obs_count == 1) {
    auto solver_args = pressio::ode::detail::tuple_take<N-1>(refs);
    auto observers   = parse_observers<Time, State>(std::get<N-1>(refs));
    return std::make_pair(std::move(solver_args), std::move(observers));
  } else {
    // all Rest... are solver args
    auto solver_args = std::move(refs);
    auto observers   = parse_observers<Time, State>();
    return std::make_pair(std::move(solver_args), std::move(observers));
  }
}

// =======================
//
// detection helpers
//
// =======================
// 4-arg explicit step
template<class S, class State, class Time, class Dt>
using explicit_step_expr_t =
  decltype(std::declval<S&>()(std::declval<State&>(),
			      std::declval<StepStartAt<Time>>(),
                              std::declval<StepCount>(),
			      std::declval<StepSize<Dt>>()));

template<class, class, class, class, class = void>
struct has_explicit_step : std::false_type {};

template<class S, class State, class Time, class Dt>
struct has_explicit_step<S, State, Time, Dt,
  std::void_t<explicit_step_expr_t<S, State, Time, Dt>>>
  : std::is_same<explicit_step_expr_t<S, State, Time, Dt>, void> {};

// 5-arg explicit step with RHS observer
template<class S, class State, class Time, class Dt, class RO>
using explicit_step_with_rhs_expr_t =
  decltype(std::declval<S&>()(std::declval<State&>(),
			      std::declval<StepStartAt<Time>>(),
                              std::declval<StepCount>(),
			      std::declval<StepSize<Dt>>(),
                              std::declval<RO&>()));

template<class, class, class, class, class, class = void>
struct has_explicit_step_with_rhs : std::false_type {};

template<class S, class State, class Time, class Dt, class RO>
struct has_explicit_step_with_rhs<S, State, Time, Dt, RO,
  std::void_t<explicit_step_with_rhs_expr_t<S, State, Time, Dt, RO>>>
  : std::is_same<explicit_step_with_rhs_expr_t<S, State, Time, Dt, RO>, void> {};


namespace detail {

// C++17 detection idiom
template<class...> using void_t = void;

template<class, template<class...> class, class...>
struct is_detected : std::false_type {};

template<template<class...> class Op, class... Args>
struct is_detected<void_t<Op<Args...>>, Op, Args...> : std::true_type {};

// Expression alias for calling impl::do_step_with_solver(..., solver, full_args...)
template<class Stepper, class State, class Time, class Dt, class Solver, class... FullArgs>
using do_step_with_solver_expr = decltype(
    std::declval<Stepper&>()(
			     std::declval<State&>(),
			     std::declval<Time>(),
			     std::declval<Dt>(),
			     std::declval<StepCount>(),
			     std::declval<Solver&>(),
			     std::declval<FullArgs&>()...
  )
);

} // namespace detail

//======================== step-size adapters =====================
namespace policy {
template<class Time>
struct fixed_dt {
  Time dt_;
  void operator()(StepCount, StepStartAt<Time>, StepSize<Time> & dt) const { dt = dt_ ; }
};
} // namespace policy

//=========================== bounds policies =====================
template<class Time, class StepSizer = policy::fixed_dt<Time>>
struct StepsFixed {
  Time t0{};
  StepCount n{};
  StepSizer step{};

  Time start_time() const { return t0; }
  bool keep_going(Time, StepCount k) const { return k <= n; }

  void set_dt(StepStartAt<Time> t, StepCount k, StepSize<Time> & dt) const {
    static_assert(std::is_void_v<decltype(step(k,t,dt))>,
                  "StepSizer must be callable with (k,t,dt) overwriting dt.");
    step(k, t, dt);
  }
};

template<class Time, class StepSizer>
struct ToTime {
  Time t0{}, tf{};
  StepSizer step{};

  Time start_time() const { return t0; }

  bool keep_going(Time t, StepCount) const {
    constexpr auto eps = std::numeric_limits<Time>::epsilon();
    if ( std::abs(t - tf) <= eps ){ return false; }
    if ( t > tf ) { return false; }
    return true;
  }

  void set_dt(StepStartAt<Time> t, StepCount k, StepSize<Time> & dt) const {
    step(k, t, dt);
  }
};


}} //end namespace pressio::ode
#endif  // PRESSIOROM_ODE_IMPL_ODE_ADVANCE_IMPL_HPP_
