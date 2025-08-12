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

#include <tuple>
#include <type_traits>
#include <utility>
#include <functional>

namespace pressio{ namespace ode{

//========================= small helpers =========================

struct state_noop_t { template<class... A> void operator()(A&&...) const noexcept {} };
struct rhs_noop_t   { template<class... A> void operator()(A&&...) const noexcept {} };

template<class T>
using obs_holder_t =
  std::conditional_t<std::is_lvalue_reference_v<T>,
    std::reference_wrapper<std::remove_reference_t<T>>,
    std::remove_reference_t<T>>;

template<class T>
constexpr auto hold(T&& v) {
  if constexpr (std::is_lvalue_reference_v<T>) return std::ref(v);
  else                                         return std::move(v);
}

//=================== observer traits + parsing ===================
template<class Time, class State, class F>
inline constexpr bool is_state_obs_v =
  StateObserver<mpl::remove_cvref_t<F>, Time, State>::value;

template<class Time, class Rhs, class F>
inline constexpr bool is_rhs_obs_v =
  RhsObserver<mpl::remove_cvref_t<F>, Time, Rhs>::value;

template<class SO, class RO>
struct parsed_obs_t { SO state; RO rhs; };

template<class Time, class State>
constexpr auto parse_observers() {
  return parsed_obs_t<state_noop_t, rhs_noop_t>{ state_noop_t{}, rhs_noop_t{} };
}

template<class Time, class State, class A1>
constexpr auto parse_observers(A1&& a1) {
  static_assert(is_state_obs_v<Time, State, A1>,
    "Single observer must be a StateObserver, placed last.");
  using SOH = obs_holder_t<A1&&>;
  return parsed_obs_t<SOH, rhs_noop_t>{ hold(std::forward<A1>(a1)), rhs_noop_t{} };
}

template<class Time, class State, class A1, class A2>
constexpr auto parse_observers(A1&& a1, A2&& a2) {
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

  constexpr bool last2_is_rhs = []() {
    if constexpr (N >= 2) {
      using L = typename std::tuple_element<N-2, Types>::type;
      return is_rhs_obs_v<Time, State, L>;
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

  void next_dt(StepStartAt<Time> t, StepCount k, StepSize<Time> & dt) const {
    static_assert(std::is_void_v<decltype(step(k,t,dt))>,
                  "StepSizer must be callable with (k,t,dt) overwriting dt.");
    step(k, t, dt);
  }
};

template<class Time, class StepSizer>
struct ToTime {
  Time t0{}, tf{};
  StepSizer step{};
  bool clip_last{true};

  Time start_time() const { return t0; }
  bool keep_going(Time t, StepCount) const {
    constexpr auto eps = std::numeric_limits<Time>::epsilon();
    if ( std::abs(t - tf) <= eps ){ return false; }
    if ( t > tf ) { return false; }
    return true;
  }

  void next_dt(StepStartAt<Time> t, StepCount k, StepSize<Time> & dt) const {
    step(k, t, dt);
    // if (clip_last && dt.get() > Time{}) {
    //   Time rem = tf - t;
    //   if (rem < dt) dt = rem;
    // }
  }
};

// Convenience factories
template<class Time>
inline StepsFixed<Time, policy::fixed_dt<Time>>
steps_fixed_dt(Time t0, StepCount n, Time dt){ return { t0, n, policy::fixed_dt<Time>{dt} }; }

template<class Time, class StepSizer>
inline StepsFixed<Time, StepSizer>
steps(Time t0, StepCount n, StepSizer sizer){ return { t0, n, std::move(sizer) }; }

template<class Time, class StepSizer>
inline ToTime<Time, StepSizer>
to_time(Time t0, Time tf, StepSizer sizer, bool clip_last=true){ return { t0, tf, std::move(sizer), clip_last }; }


//======================== advance (no solver) ====================
template<class State, class Stepper, class Policy, class... Obs>
std::enable_if_t< ExplicitStepper<Stepper>::value >
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
std::enable_if_t< !ExplicitStepper<Stepper>::value >
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
