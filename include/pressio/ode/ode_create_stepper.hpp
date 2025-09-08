/*
//@HEADER
// ************************************************************************
//
// ode_create_stepper.hpp
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

#ifndef PRESSIOROM_ODE_ODE_CREATE_STEPPER_HPP_
#define PRESSIOROM_ODE_ODE_CREATE_STEPPER_HPP_

#include "./impl/ode_user_system_wrapper.hpp"

#include "./impl/ode_explicit_stepper_without_mass_matrix.hpp"
#include "./impl/ode_explicit_stepper_with_mass_matrix.hpp"
#include "./impl/ode_implicit_discrete_residual.hpp"
#include "./impl/ode_implicit_discrete_jacobian.hpp"
#include "./impl/ode_implicit_policy_residual_jacobian_without_mass_matrix.hpp"
#include "./impl/ode_implicit_policy_residual_jacobian_with_mass_matrix.hpp"
#include "./impl/ode_implicit_stepper_standard.hpp"
#include "./impl/ode_implicit_stepper_arbitrary.hpp"
#include "./impl/ode_create_impl.hpp"

namespace pressio{ namespace ode{

namespace impl_detail {

  // unwrap unique_ptr<U,D> -> U, otherwise T unchanged
  template<class T> struct unwrap_unique_ptr { using type = T; };

  template<class U, class D>
  struct unwrap_unique_ptr<std::unique_ptr<U,D>> { using type = U; };

  template<class T>
  using unwrap_unique_ptr_t = typename unwrap_unique_ptr<T>::type;
}

/**
 * Create an explicit ODE stepper from a user-supplied system
 *
 * It accepts the system either
 *   (a) by reference (non-owning): `T&` / `const T&`
 *   (b) by transferring ownership: `std::unique_ptr<T>&&` (must be moved)
 *
 * Parameters
 * ----------
 * schemeName
 *   The explicit time-integration scheme to use.
 *   See ode_enum_and_tags.hpp for valid options to use.
 *
 * system_like
 *   The system instance or a `std::unique_ptr` to it.
 *   - Non-owning path: pass an lvalue `T&` / `const T&` (the stepper will store a reference).
 *   - Owning path:     pass `std::unique_ptr<T>` **as an rvalue** (use `std::move`);
 *                      the stepper takes ownership.
 *
 * Return
 * ------
 * An explicit stepper object.
 * Use `auto` at the call site since the type of stepper class is an impl detail.
 *
 *   MySystem sys;
 *   auto step1 = create_explicit_stepper(StepScheme::RK4, sys);   // non-owning
 *
 *   auto p = std::make_unique<MySystem>(ctor args);
 *   auto step2 = create_explicit_stepper(StepScheme::RK4, std::move(p)); // owning
 *
 * Requirements (compile-time enforced)
 * ------------------------------------
 * - `sys_type` models either `RealValuedOdeSystem` or `RealValuedOdeSystemFusingMassMatrixAndRhs`.
 * - If passing `std::unique_ptr<T>`, it must be a **non-const rvalue**
 *
 * Ownership & Lifetime
 * --------------------
 * - Non-owning path (`T&` / `const T&`): the stepper stores a reference, so
 *   the caller must ensure the system object outlives the stepper.
 * - Owning path (`std::unique_ptr<T>&&`): the stepper takes ownership of the system.
 * - Do **not** pass a temporary `T{...}` on the non-owning path (it would dangle).
 *   Wrap it in a `std::unique_ptr` and move instead.
 *
 */

template<class SystemLike>
auto create_explicit_stepper(StepScheme schemeName, SystemLike&& system_like)
{
  // Decay the parameter type: remove references and cv-qualifiers.
  // Example:
  //   SystemLike = const MySys& -> base_t = MySys
  //   SystemLike = std::unique_ptr<MySys>&& -> base_t = std::unique_ptr<MySys>
  using base_t = mpl::remove_cvref_t<SystemLike>;

  // Detect if the *decayed* type is a std::unique_ptr<...>.
  // This gates the ownership branch below (owning vs non-owning wrapper).
  constexpr bool owns = mpl::is_unique_ptr_v<base_t>;

  // Extract the underlying system type regardless of ownership.
  // unwrap_unique_ptr_t<std::unique_ptr<T>> -> T
  // unwrap_unique_ptr_t<T>                  -> T
  using sys_type = impl_detail::unwrap_unique_ptr_t<base_t>;

  // ensure the provided system exposes the expected ODE API
  static_assert(
    PRESSIO_VALUE_OF(RealValuedOdeSystem<sys_type>) ||
    PRESSIO_VALUE_OF(RealValuedOdeSystemFusingMassMatrixAndRhs<sys_type>),
    "To create an explicit stepper you must provide an ODE system that meets either "
    "RealValuedOdeSystem or RealValuedOdeSystemFusingMassMatrixAndRhs."
  );

  // get associated types from the system
  using ind_var_type = typename sys_type::independent_variable_type;
  using state_type   = typename sys_type::state_type;
  using rhs_type     = typename sys_type::rhs_type;

  // Local factory that takes a "system wrapper" (sw) and instantiates
  // the right stepper impl. W is the *decayed* wrapper type
  auto build = [&](auto&& sw){
    using W = mpl::remove_cvref_t<decltype(sw)>;

    if constexpr (PRESSIO_VALUE_OF(RealValuedOdeSystemFusingMassMatrixAndRhs<sys_type>)){
      // sys_type exposes a mass matrix type
      using mass_matrix_type = typename sys_type::mass_matrix_type;
      using impl_type = impl::ExplicitStepperWithMassMatrixImpl<
          W, state_type, ind_var_type, rhs_type, mass_matrix_type>;

      // Forward the wrapper into the impl factory
      return impl::create_explicit_stepper<impl_type>(
          schemeName, std::forward<decltype(sw)>(sw));
    } else {
      // Otherwise, use the simpler no-mass-matrix explicit stepper
      using impl_type = impl::ExplicitStepperNoMassMatrixImpl<
          W, state_type, ind_var_type, rhs_type>;

      return impl::create_explicit_stepper<impl_type>(
          schemeName, std::forward<decltype(sw)>(sw));
    }
  };

  // Ownership split
  if constexpr (owns) {
    // If the caller passed a std::unique_ptr<T>, enforce correct usage:
    // - it must be an rvalue (we are going to move-from it)
    // - it must not be const (cannot move a const unique_ptr)
    static_assert(std::is_rvalue_reference<SystemLike&&>::value &&
                  !std::is_const<std::remove_reference_t<SystemLike>>::value,
                  "Pass std::unique_ptr<T> by rvalue (use std::move).");

    // Owning wrapper: tag=1 to indicate ownership (internal convention)
    using wrapper_t = impl::SystemInternalWrapper<1, sys_type>;

    // Move the unique_ptr into the wrapper so the stepper owns the system.
    return build(wrapper_t{ std::forward<SystemLike>(system_like) });
  } else {
    // Non-owning wrapper: tag=0 to indicate reference semantics (internal convention)
    using wrapper_t = impl::SystemInternalWrapper<0, sys_type>;

    // Store by reference; the caller must ensure the system outlives the stepper.
    // Note: we intentionally do NOT forward here to avoid binding to a temporary.
    return build(wrapper_t{ system_like });
  }
}

//
// auxiliary scheme-specific functions
//
template<class ...Args>
auto create_forward_euler_stepper(Args && ...args){
  return create_explicit_stepper(StepScheme::ForwardEuler,
				 std::forward<Args>(args)...);
}

template<class ...Args>
auto create_rk4_stepper(Args && ...args){
  return create_explicit_stepper(StepScheme::RungeKutta4,
				 std::forward<Args>(args)...);
}

template<class ...Args>
auto create_ab2_stepper(Args && ...args){
  return create_explicit_stepper(StepScheme::AdamsBashforth2,
				 std::forward<Args>(args)...);
}

template<class ...Args>
auto create_ssprk3_stepper(Args && ...args){
  return create_explicit_stepper
    (StepScheme::SSPRungeKutta3, std::forward<Args>(args)...);
}





/**
 * Create an implicit ODE stepper from a user-supplied system
 *
 * It accepts the system either
 *   (a) by reference (non-owning): `T&` / `const T&`
 *   (b) by transferring ownership: `std::unique_ptr<T>&&` (must be moved)
 *
 * Parameters
 * ----------
 * schemeName
 *   The explicit time-integration scheme to use.
 *   See ode_enum_and_tags.hpp for valid options to use.
 *   Note: `CrankNicolson` is *not* supported for systems with a mass matrix
 *   (i.e., those modeling `RealValuedCompleteOdeSystem`).
 *
 * system_like
 *   The system instance or a `std::unique_ptr` to it.
 *   - Non-owning path: pass an lvalue `T&` / `const T&` (the stepper will store a reference).
 *   - Owning path:     pass `std::unique_ptr<T>` **as an rvalue** (use `std::move`);
 *                      the stepper takes ownership of the system.
 *
 * Return
 * ------
 * An implicit stepper object.
 * Use `auto` at the call site since the exact stepper type is an implementation detail.
 *
 *   MySystem sys;
 *   auto step1 = create_implicit_stepper(StepScheme::BDF2, sys);              // non-owning
 *
 *   auto p = std::make_unique<MySystem>(ctor args);
 *   auto step2 = create_implicit_stepper(StepScheme::BDF2, std::move(p));     // owning
 *
 * Requirements (compile-time enforced)
 * ------------------------------------
 * - The underlying system type models either
 *   `RealValuedOdeSystemFusingRhsAndJacobian` or `RealValuedCompleteOdeSystem`.
 * - If passing `std::unique_ptr<T>`, it must be a **non-const rvalue**.
 *
 * Ownership & Lifetime
 * --------------------
 * - Non-owning path (`T&` / `const T&`): the stepper stores a reference; the caller
 *   must ensure the system object outlives the stepper.
 * - Owning path (`std::unique_ptr<T>&&`): the stepper assumes ownership of the system.
 * - Do **not** pass a temporary `T{...}` on the non-owning path (it would dangle);
 *   wrap it in a `std::unique_ptr` and move instead.
 */

template<class SystemLike>
auto create_implicit_stepper(StepScheme schemeName, SystemLike&& system_like)
{
  // Decay the parameter type: remove references and cv-qualifiers.
  // Example:
  //   SystemLike = const MySys& -> base_t = MySys
  //   SystemLike = std::unique_ptr<MySys>&& -> base_t = std::unique_ptr<MySys>
  using base_t = mpl::remove_cvref_t<SystemLike>;

  // Detect if the *decayed* type is a std::unique_ptr<...>.
  // This gates the ownership branch below (owning vs non-owning wrapper).
  constexpr bool owns = mpl::is_unique_ptr_v<base_t>;

  // Extract the underlying system type regardless of ownership.
  // unwrap_unique_ptr_t<std::unique_ptr<T>> -> T
  // unwrap_unique_ptr_t<T>                  -> T
  using system_type = impl_detail::unwrap_unique_ptr_t<base_t>;

  // ensure the provided system exposes the expected ODE API
  constexpr bool complete_system = PRESSIO_VALUE_OF(RealValuedCompleteOdeSystem<system_type>);
  static_assert(
    PRESSIO_VALUE_OF(RealValuedOdeSystemFusingRhsAndJacobian<system_type>) ||
    complete_system,
    "To create an implicit stepper you must provide an ODE system that meets either "
    "RealValuedOdeSystemFusingRhsAndJacobian or RealValuedCompleteOdeSystem."
  );

  // Validate scheme
  assert(
    schemeName == StepScheme::BDF1 ||
    schemeName == StepScheme::BDF2 ||
    schemeName == StepScheme::CrankNicolson
  );
  if constexpr (complete_system){
    // for systems with mass matrix CN not supported yet
    assert(schemeName != StepScheme::CrankNicolson);
  }

  using ind_var_type  = typename system_type::independent_variable_type;
  using state_type    = typename system_type::state_type;
  using residual_type = typename system_type::rhs_type;
  using jacobian_type = typename system_type::jacobian_type;

  // Helper to build the policies
  auto build = [&](auto&& ws){
    using wrap_t = mpl::remove_cvref_t<decltype(ws)>;

    if constexpr (complete_system) {
      using mass_mat_type = typename system_type::mass_matrix_type;

      using policy_type = impl::ResidualJacobianWithMassMatrixStandardPolicy<
        wrap_t, ind_var_type, state_type, residual_type, jacobian_type, mass_mat_type>;

      using impl_type = impl::ImplicitStepperStandardImpl<
        ind_var_type, state_type, residual_type, jacobian_type, policy_type>;

      return impl::create_implicit_stepper_impl<impl_type>(
        schemeName, policy_type(std::forward<decltype(ws)>(ws)));
    } else {
      using policy_type = impl::ResidualJacobianStandardPolicy<
        wrap_t, ind_var_type, state_type, residual_type, jacobian_type>;

      using impl_type = impl::ImplicitStepperStandardImpl<
        ind_var_type, state_type, residual_type, jacobian_type, policy_type>;

      return impl::create_implicit_stepper_impl<impl_type>(
        schemeName, policy_type(std::forward<decltype(ws)>(ws)));
    }
  };

  if constexpr (owns) {
    // Enforce moving a unique_ptr
    static_assert(std::is_rvalue_reference<SystemLike&&>::value &&
                  !std::is_const<std::remove_reference_t<SystemLike>>::value,
                  "Pass std::unique_ptr<T> by rvalue (use std::move).");

    // Owning wrapper: tag=1 to indicate ownership (internal convention)
    using wrap_type = impl::SystemInternalWrapper<1, system_type>;
    // Move the unique_ptr into the wrapper so the stepper owns the system.
    return build(wrap_type{ std::forward<SystemLike>(system_like) });

  } else {
    // Non-owning wrapper: tag=0 to indicate reference semantics (internal convention)
    using wrap_type = impl::SystemInternalWrapper<0, system_type>;

    // Store by reference; caller must ensure lifetime outlives the stepper
    return build(wrap_type{ system_like });
  }
}

/**
 * Create an implicit ODE stepper using a *custom residual/Jacobian policy*.
 *
 * This overload lets advanced users supply their own policy that implements
 * the residual and Jacobian interfaces required by the implicit stepper.
 *
 * It accepts the policy either
 *   (a) by reference (non-owning): pass an lvalue `Policy&`
 *   (b) by value/ownership:        pass an rvalue `Policy&&` (moved into the stepper)
 *
 * Parameters
 * ----------
 * schemeName
 *   The explicit time-integration scheme to use.
 *   See ode_enum_and_tags.hpp for valid options to use.
 *
 * policy
 *   The custom residual/Jacobian policy object.
 *   - If you pass an lvalue (`policy`), the stepper will store a reference to it
 *     (non-owning). You must keep the policy alive as long as the stepper is used.
 *   - If you pass an rvalue (`std::move(policy)` or a temporary), the stepper
 *     will take it by value (moved), i.e., the stepper owns its copy.
 *
 * Return
 * ------
 * An implicit stepper object parameterized by the supplied policy.
 * Use `auto` at the call site since the concrete stepper type is an implementation detail.
 *
 * Requirements (compile-time enforced)
 * ------------------------------------
 * - `policy`’s (decayed) type models `ImplicitResidualJacobianPolicy`.
 * - The policy must provide the nested type aliases:
 *     `independent_variable_type`, `state_type`, `residual_type`, `jacobian_type`.
 *
 * Behavior
 * --------
 * - The stepper calls into the provided policy to compute the residual and Jacobian
 *   needed by the chosen implicit scheme.
 * - Scheme-specific constraints (e.g., how to handle mass matrices) are the
 *   responsibility of the custom policy.
 *
 * Ownership & Lifetime
 * --------------------
 * - Passing an lvalue (`Policy&`) makes the stepper *reference* your policy object.
 *   Ensure the policy outlives the stepper.
 * - Passing an rvalue (`Policy&&`) moves the policy into the stepper.
 *
 */

template<class ResidualJacobianPolicyType>
auto create_implicit_stepper_with_custom_policy(StepScheme schemeName,
						ResidualJacobianPolicyType && policy)
{
  using policy_type = mpl::remove_cvref_t<ResidualJacobianPolicyType>;
  static_assert(PRESSIO_VALUE_OF(ImplicitResidualJacobianPolicy<policy_type>),
  "the custom policy provided to create an implicit stepper does not meet the required concepts");

  assert(schemeName == StepScheme::BDF1 ||
	 schemeName == StepScheme::BDF2 ||
	 schemeName == StepScheme::CrankNicolson);

  using ind_var_type  = typename policy_type::independent_variable_type;
  using state_type    = typename policy_type::state_type;
  using residual_type = typename policy_type::residual_type;
  using jacobian_type = typename policy_type::jacobian_type;

  using impl_type = impl::ImplicitStepperStandardImpl<
    ind_var_type, state_type, residual_type,
    jacobian_type, ResidualJacobianPolicyType>;

  return impl::create_implicit_stepper_impl<
    impl_type>(schemeName, std::forward<ResidualJacobianPolicyType>(policy));
}

//
// auxiliary API
//
template<class ...Args>
auto create_bdf1_stepper(Args && ... args){
  return create_implicit_stepper(StepScheme::BDF1,
				 std::forward<Args>(args)...);
}

template<class ...Args>
auto create_bdf2_stepper(Args && ... args){
  return create_implicit_stepper(StepScheme::BDF2,
				 std::forward<Args>(args)...);
}

template<class ...Args>
auto create_cranknicolson_stepper(Args && ... args){
  return create_implicit_stepper(StepScheme::CrankNicolson,
				 std::forward<Args>(args)...);
}

template<class ...Args>
auto create_bdf1_stepper_with_custom_policy(Args && ... args){
  return create_implicit_stepper_with_custom_policy(StepScheme::BDF1,
						    std::forward<Args>(args)...);
}

template<class ...Args>
auto create_bdf2_stepper_with_custom_policy(Args && ... args){
  return create_implicit_stepper_with_custom_policy(StepScheme::BDF2,
						    std::forward<Args>(args)...);
}

template<class ...Args>
auto create_cranknicolson_stepper_with_custom_policy(Args && ... args){
  return create_implicit_stepper_with_custom_policy(StepScheme::CrankNicolson,
						    std::forward<Args>(args)...);
}


/**
 * Constructs a stepper for implicit time integration on a user-defined system
 * that implements Pressio’s Fully Discrete API with Jacobian support.
 *
 * Template parameters
 * - TotalNumberOfDesiredStates:
 *   Number of discrete solution states used by the stepper.
 *   Must be exactly 2 (compile-time enforced).
 *
 * Parameter
 * - odeSystem:
 *   The system to integrate. Two usage modes are supported:
 *   - Borrowed (non-owning): pass a system object or reference
 *     (e.g., create_implicit_stepper<2>(my_system);).
 *     The stepper stores a reference, the caller must keep the object alive
 *     for the entire lifetime of the stepper.
 *   - Owning: pass a std::unique_ptr<System> rvalue
 *     (e.g., create_implicit_stepper<2>(std::move(sys_ptr));).
 *     The stepper takes ownership, the pointer must be moved (not copied).
 *
 * Requirements (checked at compile time)
 * - TotalNumberOfDesiredStates == 2.
 * - Underlying system type satisfies
 *   RealValuedFullyDiscreteSystemWithJacobian<System, 2>.
 *
 * Returns
 * - An implicit stepper instance (exact type is intentionally unspecified).
 *   Use the stepper’s interface to advance the system in time.
 *
 * Examples
 * - Borrowed (non-owning):
 *     System sys{...};
 *     auto step = create_implicit_stepper<2>(sys);
 *
 * - Owning:
 *     auto sys_ptr = std::make_unique<System>(...);
 *     auto step    = create_implicit_stepper<2>(std::move(sys_ptr));
 */

template<int TotalNumberOfDesiredStates, class SystemLike>
auto create_implicit_stepper(SystemLike&& system_like)
{
  static_assert(TotalNumberOfDesiredStates == 2,
    "create_implicit_stepper currently only supports 2 total states");

  using base_t     = mpl::remove_cvref_t<SystemLike>;
  constexpr bool owns = mpl::is_unique_ptr<base_t>::value;

  using system_type = impl_detail::unwrap_unique_ptr_t<base_t>;
  static_assert(
	PRESSIO_VALUE_OF(RealValuedFullyDiscreteSystemWithJacobian<system_type,
			 TotalNumberOfDesiredStates>),
    "The ode system passed does not meet the FullyDiscrete API");

  using ind_var_type  = typename system_type::independent_variable_type;
  using state_type    = typename system_type::state_type;
  using residual_type = typename system_type::discrete_residual_type;
  using jacobian_type = typename system_type::discrete_jacobian_type;

  // Helper to build stepper without repeating the long type
  auto make_stepper = [&](auto&& sw){
    using wrap_t = mpl::remove_cvref_t<decltype(sw)>;
    using stepper_type = impl::StepperArbitrary<
        TotalNumberOfDesiredStates, wrap_t, ind_var_type,
        state_type, residual_type, jacobian_type>;
    return stepper_type(std::forward<decltype(sw)>(sw));
  };

  if constexpr (owns) {
    // Enforce rvalue for unique_ptr (must transfer ownership)
    static_assert(std::is_rvalue_reference<SystemLike&&>::value &&
                  !std::is_const<std::remove_reference_t<SystemLike>>::value,
                  "Pass std::unique_ptr<T> as an rvalue (use std::move).");

    using system_wrapper_type = impl::SystemInternalWrapper<1, system_type>;
    return make_stepper(system_wrapper_type{ std::forward<SystemLike>(system_like) });
  } else {
    // Borrowed: store a reference inside the wrapper
    using system_wrapper_type = impl::SystemInternalWrapper<0, system_type>;
    return make_stepper(system_wrapper_type{ system_like });
  }
}

}} // end namespace pressio::ode
#endif  // PRESSIOROM_ODE_ODE_CREATE_STEPPER_HPP_
