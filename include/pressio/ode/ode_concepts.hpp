/*
//@HEADER
// ************************************************************************
//
// ode_concepts.hpp
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

#ifndef PRESSIOROM_ODE_CONCEPTS_ODE_CONCEPTS_HPP_
#define PRESSIOROM_ODE_CONCEPTS_ODE_CONCEPTS_HPP_

#include "./ode_predicates.hpp"

namespace pressio{ namespace ode{

template<class T, class enable = void>
struct OdeSystem : std::false_type{};

template<class T>
struct OdeSystem<
  T,
  std::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_rhs_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::rhs_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::rhs_type >::value
    && std::is_void<
      decltype(
	       std::declval<T const>().rhs
	       (
		std::declval<typename T::state_type const&>(),
		std::declval<typename T::independent_variable_type const &>(),
		std::declval<typename T::rhs_type &>()
	       )
	   )
      >::value
   >
  > : std::true_type{};

template<class T, class enable = void>
struct OdeSystemFusingRhsAndJacobian : std::false_type{};

template<class T>
struct OdeSystemFusingRhsAndJacobian<
  T,
  std::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_rhs_typedef<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::rhs_type>::value
    && std::is_copy_constructible<typename T::jacobian_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::rhs_type >::value
    && ::pressio::ode::has_const_create_jacobian_method_return_result<
      T, typename T::jacobian_type >::value
    && std::is_void<
      decltype(
	       std::declval<T const>().rhsAndJacobian
	       (
		std::declval<typename T::state_type const&>(),
		std::declval<typename T::independent_variable_type const &>(),
                std::declval<typename T::rhs_type &>(),
		std::declval< std::optional<typename T::jacobian_type*> >()
	       )
	   )
      >::value
   >
  > : std::true_type{};


template<class T, class enable = void>
struct OdeSystemFusingMassMatrixAndRhs : std::false_type{};

template<class T>
struct OdeSystemFusingMassMatrixAndRhs<
  T,
  std::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_rhs_typedef<T>::value
    && ::pressio::has_mass_matrix_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::rhs_type>::value
    && std::is_copy_constructible<typename T::mass_matrix_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::rhs_type >::value
    && ::pressio::ode::has_const_create_mass_matrix_method_return_result<
      T, typename T::mass_matrix_type >::value
    //
    && std::is_void<
      decltype(
	       std::declval<T const>().massMatrixAndRhs
	       (
		std::declval<typename T::state_type const&>(),
		std::declval<typename T::independent_variable_type const &>(),
		std::declval<typename T::mass_matrix_type &>(),
                std::declval<typename T::rhs_type &>()
	       )
	   )
      >::value
   >
  > : std::true_type{};


template<class T, class enable = void>
struct CompleteOdeSystem : std::false_type{};

template<class T>
struct CompleteOdeSystem<
  T,
  std::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_rhs_typedef<T>::value
    && ::pressio::has_mass_matrix_typedef<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::rhs_type>::value
    && std::is_copy_constructible<typename T::mass_matrix_type>::value
    && std::is_copy_constructible<typename T::jacobian_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::rhs_type >::value
    && ::pressio::ode::has_const_create_mass_matrix_method_return_result<
      T, typename T::mass_matrix_type >::value
    && ::pressio::ode::has_const_create_jacobian_method_return_result<
      T, typename T::jacobian_type >::value
    //
    && std::is_void<
      decltype(
	       std::declval<T const>().massMatrixAndRhsAndJacobian
	       (
		std::declval<typename T::state_type const&>(),
		std::declval<typename T::independent_variable_type const &>(),
                std::declval<typename T::mass_matrix_type &>(),
                std::declval<typename T::rhs_type &>(),
		std::declval< std::optional<typename T::jacobian_type*> >()
	       )
	   )
      >::value
   >
  > : std::true_type{};

template<class T, int NumStates,class enable = void>
struct FullyDiscreteSystemWithJacobian : std::false_type{};

template<class T, int NumStates>
struct FullyDiscreteSystemWithJacobian<
  T, NumStates,
  std::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_discrete_residual_typedef<T>::value
    && ::pressio::has_discrete_jacobian_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::discrete_residual_type>::value
    && std::is_copy_constructible<typename T::discrete_jacobian_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_discrete_residual_method_return_result<
      T, typename T::discrete_residual_type>::value
    && ::pressio::ode::has_const_create_discrete_jacobian_method_return_result<
      T, typename T::discrete_jacobian_type>::value
    //
    && ::pressio::ode::has_const_discrete_residual_jacobian_method<
      T, NumStates,
      typename ::pressio::ode::StepCount::value_type,
      typename T::independent_variable_type,
      typename T::state_type,
      typename T::discrete_residual_type,
      typename T::discrete_jacobian_type
      >::value
    >
  > : std::true_type{};


//
// refine for real-valued case
//
template<class T, class enable = void>
struct RealValuedOdeSystem : std::false_type{};

template<class T>
struct RealValuedOdeSystem<
  T, std::enable_if_t<
       OdeSystem<T>::value
       && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
       && std::is_floating_point< scalar_trait_t<typename T::rhs_type> >::value
       && std::is_convertible<
	 typename T::independent_variable_type,
	 scalar_trait_t<typename T::state_type> >::value
  >
  > : std::true_type{};

template<class T, class enable = void>
struct RealValuedOdeSystemFusingRhsAndJacobian : std::false_type{};

template<class T>
struct RealValuedOdeSystemFusingRhsAndJacobian<
  T,
  std::enable_if_t<
    OdeSystemFusingRhsAndJacobian<T>::value
  && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
  && std::is_floating_point< scalar_trait_t<typename T::rhs_type> >::value
  && std::is_floating_point< scalar_trait_t<typename T::jacobian_type> >::value
  && std::is_convertible<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >::value
  > > : std::true_type{};


template<class T, class enable = void>
struct RealValuedOdeSystemFusingMassMatrixAndRhs : std::false_type{};

template<class T>
struct RealValuedOdeSystemFusingMassMatrixAndRhs<
  T, std::enable_if_t<
    OdeSystemFusingMassMatrixAndRhs<T>::value
  && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
  && std::is_floating_point< scalar_trait_t<typename T::rhs_type> >::value
  && std::is_floating_point< scalar_trait_t<typename T::mass_matrix_type> >::value
  && std::is_convertible<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >::value
  >
  > : std::true_type{};


template<class T, class enable = void>
struct RealValuedCompleteOdeSystem : std::false_type{};

template<class T>
struct RealValuedCompleteOdeSystem<
  T, std::enable_if_t<
       CompleteOdeSystem<T>::value
       && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
       && std::is_floating_point< scalar_trait_t<typename T::rhs_type> >::value
       && std::is_floating_point< scalar_trait_t<typename T::mass_matrix_type> >::value
       && std::is_floating_point< scalar_trait_t<typename T::jacobian_type> >::value
       && std::is_convertible<
	 typename T::independent_variable_type,
	 scalar_trait_t<typename T::state_type> >::value
  >
  > : std::true_type{};

template <class T, int NumStates, class = void>
struct RealValuedFullyDiscreteSystemWithJacobian : std::false_type{};

template <class T, int NumStates>
struct RealValuedFullyDiscreteSystemWithJacobian<
  T, NumStates,
  std::enable_if_t<
    FullyDiscreteSystemWithJacobian<T, NumStates>::value
    && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
    && std::is_floating_point< scalar_trait_t<typename T::discrete_residual_type> >::value
    && std::is_floating_point< scalar_trait_t<typename T::discrete_jacobian_type> >::value
    && std::is_convertible<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >::value
  >
  > : std::true_type{};


//
// policy
//
template<class T, class = void>
struct ImplicitResidualJacobianPolicy : std::false_type{};

template<class T>
struct ImplicitResidualJacobianPolicy<
  T,
  std::enable_if_t<
    ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_residual_typedef<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    //
    && ::pressio::ops::is_known_data_type<typename T::state_type>::value
    && ::pressio::ops::is_known_data_type<typename T::residual_type>::value
    && ::pressio::ops::is_known_data_type<typename T::jacobian_type>::value
    && all_have_traits_and_same_scalar<
      typename T::state_type,
      typename T::residual_type,
      typename T::jacobian_type>::value
    && std::is_convertible<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type>>::value
    //
    // create methods
    //
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type>::value
    && std::is_same<
      typename T::residual_type,
      decltype(std::declval<T const>().createResidual())
      >::value
    && std::is_same<
      typename T::jacobian_type,
      decltype(std::declval<T const>().createJacobian())
      >::value

    && std::is_void<
      decltype
      (
       std::declval<T const>()
       (
	std::declval<StepScheme const &>(),
	std::declval<typename T::state_type const &>(),
	std::declval<ImplicitStencilStatesDynamicContainer<typename T::state_type> const & >(),
	std::declval<ImplicitStencilRightHandSideDynamicContainer<typename T::residual_type> & >(),
	std::declval< ::pressio::ode::StepEndAt<typename T::independent_variable_type> >(),
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepSize<typename T::independent_variable_type> >(),
	std::declval<typename T::residual_type &>(),
	std::declval< std::optional<typename T::jacobian_type*> >()
	)
       )
      >::value
    >
  > : std::true_type{};


/*
  This comment documents the pieces defined just below:

    1) `stepper_call_t<T>` — an alias for the *deduced return type* of calling
       `T::operator()` with the expected arguments.
    2) `StepperWithoutSolver<T>` — a boolean type trait that evaluates to
       `true` iff that call is well-formed and returns **exactly `void`**.

  ---------------------------------------------------------------------------
  What does the trait check?
  ---------------------------------------------------------------------------
  Given a `T` modeling a stepper, we expect:

    - `T::independent_variable_type`   // e.g., double (the "time" type)
    - `T::state_type`                  // the state container to be advanced
    - A callable operator with the shape:

        void operator()( state_type&,
                         StepStartAt<independent_variable_type>,
                         StepCount,
                         StepSize<independent_variable_type> );

  If such a call is viable and its deduced return type is **exactly `void`**,
  then `StepperWithoutSolver<T>::value == true`; otherwise `false`.

  ---------------------------------------------------------------------------
  How does it work?
  ---------------------------------------------------------------------------
  1) `stepper_call_t<T>` forms the type of the expression

        std::declval<T&>()(
          std::declval<typename T::state_type&>(),
          std::declval<StepStartAt<typename T::independent_variable_type>>(),
          std::declval<StepCount>(),
          std::declval<StepSize<typename T::independent_variable_type>>()
        )

     using `decltype(...)`.  `std::declval<X>()` fabricates a value of type `X`
     in an unevaluated context, letting us ask “what would the call return?”
     If the call is ill-formed (wrong/missing types), `stepper_call_t<T>` itself
     fails to form and is SFINAE’d out.

  2) The partial specialization of `StepperWithoutSolver<T, ...>` is gated by
     `std::void_t<...>`:

        std::void_t<
          typename T::independent_variable_type,
          typename T::state_type,
          stepper_call_t<T>
        >

     This specialization only participates if the nested typedefs exist *and*
     the call expression type could be formed.

  3) Inside the participating specialization we inherit from
     `std::is_same<void, stepper_call_t<T>>`, which makes the trait evaluate to
     `std::true_type` only when the deduced return type is **exactly `void`**.

  ---------------------------------------------------------------------------
  Enforced:
    - `T` provides `state_type` and `independent_variable_type`.
    - The call with the shown argument types is viable.
    - The return type is **exactly `void`**

  ---------------------------------------------------------------------------
  Usage examples
  ---------------------------------------------------------------------------
  struct Good {
    using independent_variable_type = double;
    using state_type = State;
    void operator()(state_type&, StepStartAt<double>, StepCount, StepSize<double>) {}
  };
  static_assert(StepperWithoutSolver<Good>::value, "Good stepper");

  struct BadReturn {
    using independent_variable_type = double;
    using state_type = State;
    int operator()(state_type&, StepStartAt<double>, StepCount, StepSize<double>) { return 0; }
  };
  static_assert(!StepperWithoutSolver<BadReturn>::value, "must return void");
*/

template<class T>
using stepper_call_t = decltype(
  std::declval<T&>()(
    std::declval<typename T::state_type &>(),
    std::declval<StepStartAt<typename T::independent_variable_type>>(),
    std::declval<StepCount>(),
    std::declval<StepSize<typename T::independent_variable_type>>()
  )
);

template <class T, class = void>
struct StepperWithoutSolver : std::false_type{};

template <class T>
struct StepperWithoutSolver<
  T,
  std::void_t<
    typename T::independent_variable_type,
    typename T::state_type,
    stepper_call_t<T>
    >
>
: std::is_same<void, stepper_call_t<T>> {};

// -------------------------------------------------------------
// -------------------------------------------------------------

/*
  Trait to detect a **const** state observer callable.

  Succeeds iff a type `T` has a **const** call operator with the shape:

    void operator()( ::pressio::ode::StepCount,
                     IndVarType,
                     StateType const& ) const;

  More precisely:
   - We form the type of the expression
       std::declval<T const&>()(StepCount{}, IndVarType{}, StateType const&{})
     via `decltype(...)`. If that call is ill-formed, the specialization
     is discarded (SFINAE) and the trait is `false`.
   - If the call is well-formed, we additionally require the deduced
     return type to be **exactly `void`** (not just convertible to void).

  Usage:
    static_assert(StateObserver<MyObs, double, MyState>::value, "Must be a const observer");
*/

template<class T, class IndVarType, class StateType>
using state_observer_call_c_t = decltype(
  std::declval<T const&>()(
    std::declval<StepCount>(),
    std::declval<IndVarType>(),
    std::declval<StateType const&>()
  )
);

template<class T, class IndVarType, class StateType, class = void>
struct StateObserver : std::false_type {};

template<class T, class IndVarType, class StateType>
struct StateObserver<
  T, IndVarType, StateType,
  std::void_t< state_observer_call_c_t<T, IndVarType, StateType> >
>
: std::is_same<void, state_observer_call_c_t<T, IndVarType, StateType>> {};




// -------------------------------------------------------------
// -------------------------------------------------------------

// Detect a **const** RHS observer with exact-void return:
//   void operator()(StepCount, IntermediateStepCount, IndVarType, RhsType const&) const;

template<class T, class IndVarType, class RhsType>
using rhs_observer_call_c_t = decltype(
  std::declval<T const&>()(
    std::declval<StepCount>(),
    std::declval<IntermediateStepCount>(),
    std::declval<IndVarType>(),
    std::declval<RhsType const&>()
  )
);

// Primary: default to false
template<class T, class IndVarType, class RhsType, class = void>
struct RhsObserver : std::false_type {};

// Specialization: only if the const-call is well-formed; require return type == void
template<class T, class IndVarType, class RhsType>
struct RhsObserver<
  T, IndVarType, RhsType,
  std::void_t< rhs_observer_call_c_t<T, IndVarType, RhsType> >
>
: std::is_same<void, rhs_observer_call_c_t<T, IndVarType, RhsType>> {};



// -------------------------------------------------------------
// -------------------------------------------------------------

/*
  Trait: StepSizePolicy<T, IndVarType>

  Goal
  ----
  Detect at compile time whether a type `T` can act as a *step size policy*,
  i.e., it is a callable object with a **const** call operator of the form:

      void operator()( ::pressio::ode::StepCount,
                       ::pressio::ode::StepStartAt<IndVarType>,
                       ::pressio::ode::StepSize<IndVarType>& ) const;

  Semantics
  ---------
  - The third argument is a non-const reference to StepSize<IndVarType>, so the
    policy can write/adjust the step size.
  - The trait requires the return type to be **exactly `void`**.
  - If the call expression is ill-formed (wrong/missing arguments), the trait
    evaluates to false (no hard compile error by itself).

  How it works (brief)
  --------------------
  1) `step_size_policy_call_c_t<...>` forms the type of the expression
         std::declval<T const&>()(StepCount, StepStartAt<IndVarType>, StepSize<IndVarType>&)
     using `decltype(...)`. If this fails, SFINAE removes the specialization.
  2) The partial specialization of `StepSizePolicy` only participates if that
     call type exists. It then inherits from `std::is_same<void, ...>` to
     enforce an exact `void` return.

  Notes
  -----
  - This version expects a **const** call operator.
*/

template<class T, class IndVarType>
using step_size_policy_call_c_t = decltype(
  std::declval<T const&>()(
    std::declval<StepCount>(),
    std::declval<StepStartAt<IndVarType>>(),
    std::declval<StepSize<IndVarType>&>()
  )
);

template <class T, class IndVarType, class Enable = void>
struct StepSizePolicy : std::false_type {};

template <class T, class IndVarType>
struct StepSizePolicy<
  T, IndVarType,
  std::void_t< step_size_policy_call_c_t<T, IndVarType> >
>
: std::is_same<void, step_size_policy_call_c_t<T, IndVarType>> {};


// -------------------------------------------------------------
// -------------------------------------------------------------

/*
  Trait: StepSizePolicyWithReductionScheme<T, IndVarType>

  Detect at compile time whether a type `T` can act as a *step size policy with
  reduction scheme*, i.e., it has a **const** call operator with the shape:

      void operator()( ::pressio::ode::StepCount,
                       ::pressio::ode::StepStartAt<IndVarType>,
                       ::pressio::ode::StepSize<IndVarType>&,
                       ::pressio::ode::StepSizeMinAllowedValue<IndVarType>&,
                       ::pressio::ode::StepSizeScalingFactor<IndVarType>& ) const;

  Semantics
  ---------
  - The last three parameters are non-const references the policy can modify:
      * StepSize<IndVarType>&                — current step size
      * StepSizeMinAllowedValue<IndVarType>& — minimum allowed step size
      * StepSizeScalingFactor<IndVarType>&   — multiplicative scaling factor
  - The trait requires the return type to be **exactly `void`**.

  Notes
  -----
  - This detects a **const** call operator. If you also want to accept non-const,
    add a parallel detection with `T&` and OR the results.
*/

template<class T, class IndVarType>
using step_size_policy_with_reduction_call_c_t = decltype(
  std::declval<T const&>()(
    std::declval<StepCount>(),
    std::declval<StepStartAt<IndVarType>>(),
    std::declval<StepSize<IndVarType>&>(),
    std::declval<StepSizeMinAllowedValue<IndVarType>&>(),
    std::declval<StepSizeScalingFactor<IndVarType>&>()
  )
);

template <class T, class IndVarType, class Enable = void>
struct StepSizePolicyWithReductionScheme : std::false_type {};

template <class T, class IndVarType>
struct StepSizePolicyWithReductionScheme<
  T, IndVarType,
  std::void_t< step_size_policy_with_reduction_call_c_t<T, IndVarType> >
>
: std::is_same<void, step_size_policy_with_reduction_call_c_t<T, IndVarType>> {};

}}
#endif  // PRESSIOROM_ODE_CONCEPTS_ODE_SYSTEM_HPP_
