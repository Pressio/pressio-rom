/*
//@HEADER
// ************************************************************************
//
// ode_create_explicit_stepper.hpp
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

#ifndef PRESSIOROM_ODE_ODE_CREATE_EXPLICIT_STEPPER_HPP_
#define PRESSIOROM_ODE_ODE_CREATE_EXPLICIT_STEPPER_HPP_

#include "./impl/ode_user_system_wrapper.hpp"
#include "./impl/ode_explicit_stepper_without_mass_matrix.hpp"
#include "./impl/ode_explicit_stepper_with_mass_matrix.hpp"
#include "./impl/ode_explicit_create_impl.hpp"

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
 * Template Parameters
 * -------------------
 * SystemLike
 *   A forwarding type that is either `T`, `T&`, `const T&`, `std::unique_ptr<T>`, or
 *   `std::unique_ptr<T>&&`. CV-qualifiers and references are removed for trait checks.
 *
 * Parameters
 * ----------
 * schemeName
 *   The explicit time-integration scheme to use (e.g., Forward Euler, RK3, RK4).
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
 * Use `auto` at the call site:
 *
 *   MySystem sys;
 *   auto step1 = create_explicit_stepper(StepScheme::RK4, sys);   // non-owning
 *
 *   auto p = std::make_unique<MySystem>(/* ctor args * /);
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

  // ensure the provided system exposes the expected ODE API.
  static_assert(
    RealValuedOdeSystem<sys_type>::value ||
    RealValuedOdeSystemFusingMassMatrixAndRhs<sys_type>::value,
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

    if constexpr (RealValuedOdeSystemFusingMassMatrixAndRhs<sys_type>::value){
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

}} // end namespace pressio::ode
#endif  // PRESSIOROM_ODE_ODE_CREATE_EXPLICIT_STEPPER_HPP_
