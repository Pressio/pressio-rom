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

// ---------- small helpers ----------

template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

struct state_noop_t { template<class... A> void operator()(A&&...) const noexcept {} };
struct rhs_noop_t   { template<class... A> void operator()(A&&...) const noexcept {} };

// keep references for lvalues; move for rvalues
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

// ---------- dt helper (C++17) ----------

template<class F, class Time>
struct is_dt_like : std::bool_constant<
  std::is_arithmetic_v<std::decay_t<F>> ||
  std::is_convertible_v<std::decay_t<F>, Time>
>{};

// ---------- observer traits (use your existing StateObserver/RhsObserver) ----------

template<class Time, class State, class F>
inline constexpr bool is_state_obs_v =
  StateObserver<mpl::remove_cvref_t<F>, Time, State>::value;

template<class Time, class Rhs, class F>
inline constexpr bool is_rhs_obs_v =
  RhsObserver<mpl::remove_cvref_t<F>, Time, Rhs>::value;

// ---------- parsed observers pack ----------

template<class SO, class RO, std::size_t Consumed>
struct parsed_obs_t {
  SO state;
  RO rhs;
  static constexpr std::size_t consumed = Consumed;
};

template<class Time, class State>
constexpr auto parse_observers() {
  return parsed_obs_t<state_noop_t, rhs_noop_t, 0>{ state_noop_t{}, rhs_noop_t{} };
}

template<class Time, class State, class A1>
constexpr auto parse_observers(A1&& a1) {
  static_assert(is_state_obs_v<Time, State, A1>,
    "Single observer must be a StateObserver, placed last.");
  using SOH = obs_holder_t<A1&&>;
  return parsed_obs_t<SOH, rhs_noop_t, 1>{ hold(std::forward<A1>(a1)), rhs_noop_t{} };
}

template<class Time, class State, class A1, class A2>
constexpr auto parse_observers(A1&& a1, A2&& a2) {
  static_assert(is_state_obs_v<Time, State, A1> && is_rhs_obs_v<Time, State, A2>,
    "Two observers must be (StateObserver, RhsObserver), placed last.");
  using SOH = obs_holder_t<A1&&>;
  using ROH = obs_holder_t<A2&&>;
  return parsed_obs_t<SOH, ROH, 2>{ hold(std::forward<A1>(a1)), hold(std::forward<A2>(a2)) };
}

// ---------- split (solver args..., observers) from the tail ----------

namespace detail {

// true if any of the first K tuple elements is an observer
template<class Time, class State, class Tup, std::size_t... I>
constexpr bool any_observer_in_prefix(std::index_sequence<I...>) {
  return ((is_state_obs_v<Time, State, std::remove_reference_t<std::tuple_element_t<I, Tup>>> ||
           is_rhs_obs_v  <Time, State, std::remove_reference_t<std::tuple_element_t<I, Tup>>>) || ...);
}

template<std::size_t K, class Tup, std::size_t... I>
constexpr auto tuple_take_impl(Tup&& tup, std::index_sequence<I...>) {
  return std::forward_as_tuple(std::get<I>(std::forward<Tup>(tup))...);
}
template<std::size_t K, class Tup>
constexpr auto tuple_take(Tup&& tup) {
  if constexpr (K == 0) return std::tuple<>{};
  else                  return tuple_take_impl<K>(std::forward<Tup>(tup), std::make_index_sequence<K>{});
}

} // namespace detail

template<class Time, class State, class Solver, class... Rest>
constexpr auto split_solver_and_observers(Solver&, Rest&&... rest)
{
  using Tup = std::tuple<Rest&&...>;
  Tup tup(std::forward<Rest>(rest)...);
  constexpr std::size_t N = sizeof...(Rest);

  // decide trailing observers: check last two, then last one
  constexpr bool last_is_state = (N >= 1) ?
    is_state_obs_v<Time, State, std::remove_reference_t<std::tuple_element_t<N-1, Tup>>> : false;
  constexpr bool last_is_rhs   = (N >= 1) ?
    is_rhs_obs_v  <Time, State, std::remove_reference_t<std::tuple_element_t<N-1, Tup>>> : false;

  constexpr bool last2_is_state = (N >= 2) ?
    is_state_obs_v<Time, State, std::remove_reference_t<std::tuple_element_t<N-2, Tup>>> : false;
  constexpr bool last2_is_rhs   = (N >= 2) ?
    is_rhs_obs_v  <Time, State, std::remove_reference_t<std::tuple_element_t<N-2, Tup>>> : false;

  static_assert(!(last_is_rhs && !(N >= 2 && last2_is_state)),
    "RhsObserver cannot appear without a preceding StateObserver.");

  constexpr std::size_t obs_count =
    (N >= 2 && last2_is_state && last_is_rhs) ? 2u :
    (N >= 1 && last_is_state && !last_is_rhs) ? 1u : 0u;

  // ensure observers are trailing only
  if constexpr (obs_count <= N) {
    constexpr std::size_t K = N - obs_count;
    static_assert(!detail::any_observer_in_prefix<Time, State, Tup>(std::make_index_sequence<K>{}),
      "Observers must be the last arguments: none, (StateObserver), or (StateObserver,RhsObserver).");
  }

  // solver args tuple (prefix)
  auto solver_args_tuple = [&](){
    if constexpr (obs_count == 2)      return detail::tuple_take<N-2>(tup);
    else if constexpr (obs_count == 1) return detail::tuple_take<N-1>(tup);
    else                                return tup;
  }();

  // parsed observers
  auto observers = [&](){
    if constexpr (obs_count == 2)
      return parse_observers<Time, State>(std::get<N-2>(tup), std::get<N-1>(tup));
    else if constexpr (obs_count == 1)
      return parse_observers<Time, State>(std::get<N-1>(tup));
    else
      return parse_observers<Time, State>();
  }();

  return std::make_pair(std::move(solver_args_tuple), std::move(observers));
}

// ---------- advance() #1: no solver (only observers at the end) ----------

template<
  class StepperType,
  class StateType,
  class IndVarType,
  class FirstArg,
  class... Rest>
std::enable_if_t< ExplicitStepper<StepperType>::value >
advance(StepperType & stepper,
        StateType   & state,
        const IndVarType & startVal,
        FirstArg && first,                       // dt or StepSizePolicy
        std::variant<StepCount, IndVarType> var, // steps or final time
        Rest&&... rest)                          // observers only
{
  using Time = IndVarType;

  static_assert(sizeof...(Rest) <= 2,
  "Only observers are allowed here: none, (StateObserver), or (StateObserver,RhsObserver).");

  auto obs = parse_observers<Time, StateType>(std::forward<Rest>(rest)...);
  auto&& state_obs = obs.state;
  auto&& rhs_obs   = obs.rhs;

  Time t = startVal;
  state_obs(StepCount{0}, t, state);

  StepCount k = StepCount{::pressio::ode::first_step_value};

  // auto do_one_with_dt = [&](Time dt){
  //   stepper(state, t, dt, k);
  //   t += dt; ++k;
  //   state_obs(k, t, state);
  // };

  auto run_fixed_dt = [&](Time dt){
    std::visit(overloaded{
      [&](StepCount n){
	for (StepCount i=0; i<n; ++i) do_one_with_dt(dt);
      },
      [&](Time tf){
	while (t < tf) do_one_with_dt(dt);
      }
    }, var);
  };

  // auto run_policy = [&](auto&& policy){
  //   std::visit(overloaded{
  //     [&](StepCount n){
  //       for (StepCount i=0; i<n; ++i){
  //         Time pdt = policy.next_dt(state, t, k);
  //         do_one_with_dt(pdt);
  //       }
  //     },
  //     [&](Time tf){
  //       while (t < tf){
  //         Time pdt = policy.next_dt(state, t, k);
  //         do_one_with_dt(pdt);
  //       }
  //     }
  //   }, var);
  // };

  if constexpr (is_dt_like<FirstArg, Time>::value) {
    run_fixed_dt(static_cast<Time>(first));
  }
  // else if constexpr (StepSizePolicy<std::decay_t<FirstArg>, StateType, Time>::value) { // <—
  //   run_policy(std::forward<FirstArg>(first));
  // } else {
  //   static_assert([]{return false;}(),
  //     "The 5th argument must be either a fixed dt (convertible to Time) "
  //     "or a StepSizePolicy with next_dt(state,t,k).");
  // }
}

template<
  class StepperType,
  class StateType,
  class IndVarType,
  class FirstArg,
  class... Rest>
inline void
advance(StepperType & stepper,
        StateType   & state,
        const IndVarType & startVal,
        FirstArg && first,
        StepCount sc,
        Rest&&... rest)
{
  advance(stepper, state, startVal, std::forward<FirstArg>(first),
	  std::variant<StepCount,IndVarType>{sc},
	  std::forward<Rest>(rest)...);
}

template<
  class StepperType,
  class StateType,
  class IndVarType,
  class FirstArg,
  class... Rest>
inline void
advance(StepperType & stepper,
        StateType   & state,
        const IndVarType & startVal,
        FirstArg && first,
        IndVarType tf,
        Rest&&... rest)
{
  advance(stepper, state, startVal, std::forward<FirstArg>(first),
	  std::variant<StepCount,IndVarType>{tf},
	  std::forward<Rest>(rest)...);
}

// ---------- advance() #2: with solver (solver args..., [SO], [SO,RO]) ----------

template<
  class StepperType, class StateType, class IndVarType,
  class FirstArg, class SolverType, class... Rest>
void advance(StepperType & stepper,
             StateType   & state,
             const IndVarType & startVal,
             std::variant<StepCount, IndVarType> var, // steps or final time
             FirstArg && first,                       // dt or StepSizePolicy
             SolverType & solver,
             Rest&&... rest)                          // solver args..., observers
{
  // using Time = IndVarType;

  // auto split = split_solver_and_observers<Time, StateType>(solver, std::forward<Rest>(rest)...);
  // auto& solver_args = split.first;
  // auto  obs         = std::move(split.second);

  // auto&& state_obs = obs.state;
  // auto&& rhs_obs   = obs.rhs;

  // Time t = startVal;
  // StepCount k = 1;

  // auto do_one_with_dt = [&](Time dt){
  //   auto call = [&](auto&... args){
  //     // default hook; you implement this in impl::
  //     impl::do_step_with_solver(stepper, state, t, dt, k, solver, args...);
  //   };
  //   std::apply(call, solver_args);

  //   state_obs(k, t, state);
  //   if constexpr (requires{ stepper.last_rhs(); })
  //     rhs_obs(k, t, state, stepper.last_rhs());
  //   t += dt; ++k;
  // };

  // auto run_fixed_dt = [&](Time dt){
  //   std::visit(overloaded{
  //     [&](StepCount n){ for (StepCount i=0; i<n; ++i) do_one_with_dt(dt); },
  //     [&](Time tf){ while (t < tf) do_one_with_dt(dt); }
  //   }, var);
  // };

  // auto run_policy = [&](auto&& policy){
  //   std::visit(overloaded{
  //     [&](StepCount n){
  //       for (StepCount i=0; i<n; ++i){
  //         Time pdt = policy.next_dt(state, t, k);
  //         do_one_with_dt(pdt);
  //       }
  //     },
  //     [&](Time tf){
  //       while (t < tf){
  //         Time pdt = policy.next_dt(state, t, k);
  //         do_one_with_dt(pdt);
  //       }
  //     }
  //   }, var);
  // };

  // if constexpr (is_dt_like<FirstArg, Time>::value) {
  //   run_fixed_dt(static_cast<Time>(first));
  // } else if constexpr (StepSizePolicy<std::decay_t<FirstArg>, StateType, Time>::value) { // <—
  //   run_policy(std::forward<FirstArg>(first));
  // } else {
  //   static_assert([]{return false;}(),
  //     "The 5th argument must be either a fixed dt (convertible to Time) "
  //     "or a StepSizePolicy with next_dt(state,t,k).");
  // }
}






// template<
//   class StepperType,
//   class StateType,
//   class IndVarType,
//   class FirstArg,        // either const dt (IndVarType) or a StepSizePolicy
//   class... Rest          // 0, 1, or 2: [StateObserver], [StateObserver, RhsObserver]
// >
// std::enable_if_t< ExplicitStepper<StepperType>::value >
// advance_n_steps(StepperType & stepper,
//                 StateType & state,
//                 const IndVarType & startVal,
//                 FirstArg && first,
//                 StepCount numSteps,
//                 Rest&&... rest)
// {
//   static_assert(sizeof...(Rest) <= 2,
//     "advance_n_steps: only (no observers), (StateObserver), or (StateObserver, RhsObserver) are allowed.");

//   impl::mandate_on_ind_var_and_state_types(stepper, state, startVal);

//   using state_noop_t = impl::NoOpStateObserver<IndVarType, StateType>;
//   using rhs_noop_t   = impl::NoOpRhsObserver<IndVarType, StateType>;

//   // Parse optional observers
//   auto make_observers = [] (auto&&... tail)
//   {
//     // if the pack is empty, it means we don't want any observer
//     if constexpr (sizeof...(tail) == 0) {
//       return std::pair<state_noop_t, rhs_noop_t>{};
//     }
//     // if there is only one, then it must be a state observer
//     else if constexpr (sizeof...(tail) == 1) {
//       auto && s = (std::forward<decltype(tail)>(tail), ...);
//       static_assert(StateObserver<mpl::remove_cvref_t<decltype(s)>, IndVarType, StateType>::value,
//         "advance_n_steps: single trailing argument must be a valid StateObserver.");

//       return std::pair<decltype(s), rhs_noop_t>{ std::forward<decltype(s)>(s), rhs_noop_t{} };
//     }
//     else { // exactly 2
//       auto tup = std::forward_as_tuple(std::forward<decltype(tail)>(tail)...);
//       auto && s = std::get<0>(tup);
//       auto && r = std::get<1>(tup);
//       static_assert(StateObserver<mpl::remove_cvref_t<decltype(s)>, IndVarType, StateType>::value,
//         "advance_n_steps: first trailing argument must be a valid StateObserver.");

//       static_assert(RhsObserver<mpl::remove_cvref_t<decltype(r)>, IndVarType, StateType>::value,
//         "advance_n_steps: second trailing argument must be a valid RhsObserver.");
//       return std::pair<decltype(s), decltype(r)>{ std::forward<decltype(s)>(s),
//                                                   std::forward<decltype(r)>(r) };
//     }
//   };
//   auto [stateObs, rhsObs] = make_observers(std::forward<Rest>(rest)...);

//   // Dispatch on policy vs constant dt
//   if constexpr (StepSizePolicy<FirstArg&&, IndVarType>::value) {
//     impl::advance_n_steps_with_dt_policy(stepper,
//                                          numSteps,
//                                          startVal,
//                                          state,
//                                          std::forward<FirstArg>(first),
//                                          std::forward<decltype(stateObs)>(stateObs),
//                                          std::forward<decltype(rhsObs)>(rhsObs));
//   } else {
//     static_assert(std::is_convertible_v<std::decay_t<FirstArg>, IndVarType>,
//       "advance_n_steps: when not passing a StepSizePolicy, `first` must be convertible to IndVarType (step size).");
//     const IndVarType & stepSize = first;
//     impl::advance_n_steps_with_dt_policy(stepper,
// 					 numSteps,
// 					 startVal,
// 					 state,
// 					 [sz = stepSize](StepCount /*currStep*/,
// 							 StepStartAt<IndVarType> /*currTime*/,
// 							 StepSize<IndVarType> & dt){
// 					   dt = sz;
// 					 },
// 					 std::forward<decltype(stateObs)>(stateObs),
// 					 std::forward<decltype(rhsObs)>(rhsObs));
//   }
// }


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
