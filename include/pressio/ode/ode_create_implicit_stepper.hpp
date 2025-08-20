/*
//@HEADER
// ************************************************************************
//
// ode_create_implicit_stepper.hpp
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

#ifndef PRESSIOROM_ODE_ODE_CREATE_IMPLICIT_STEPPER_HPP_
#define PRESSIOROM_ODE_ODE_CREATE_IMPLICIT_STEPPER_HPP_


#include "./impl/ode_implicit_discrete_residual.hpp"
#include "./impl/ode_implicit_discrete_jacobian.hpp"
#include "./impl/ode_user_system_wrapper.hpp"
#include "./impl/ode_implicit_policy_residual_jacobian_without_mass_matrix.hpp"
#include "./impl/ode_implicit_policy_residual_jacobian_with_mass_matrix.hpp"
#include "./impl/ode_implicit_stepper_standard.hpp"
#include "./impl/ode_implicit_stepper_arbitrary.hpp"
#include "./impl/ode_implicit_create_impl.hpp"

namespace pressio{ namespace ode{

// template<class SystemType>
// auto create_implicit_stepper(StepScheme schemeName,
// 			     SystemType const & odeSystem)
// {
//   using system_type = mpl::remove_cvref_t<SystemType>;
//   static_assert(RealValuedOdeSystemFusingRhsAndJacobian<system_type>::value ||
// 		RealValuedCompleteOdeSystem<system_type>::value,
// 		"To create an implicit stepper you must provide an ODE system that meets either the RealValuedOdeSystemFusingRhsAndJacobian or the RealValuedCompleteOdeSystem concept.");
//   static constexpr bool complete_system = RealValuedCompleteOdeSystem<system_type>::value;

//   assert(
// 	 schemeName == StepScheme::BDF1 ||
// 	 schemeName == StepScheme::BDF2 ||
// 	 schemeName == StepScheme::CrankNicolson
// 	 );
//   if constexpr(complete_system){
//     // for systems with mass matrix CN not supported yet
//     assert(schemeName != StepScheme::CrankNicolson);
//   }

//   using ind_var_type  = typename system_type::independent_variable_type;
//   using state_type    = typename system_type::state_type;
//   using residual_type = typename system_type::rhs_type;
//   using jacobian_type = typename system_type::jacobian_type;

//   using wrap_type = impl::SystemInternalWrapper<0, system_type>;
//   wrap_type ws(odeSystem);

//   if constexpr (complete_system){
//     using mass_mat_type = typename system_type::mass_matrix_type;

//     using policy_type = impl::ResidualJacobianWithMassMatrixStandardPolicy<
//       wrap_type, ind_var_type, state_type,
//       residual_type, jacobian_type, mass_mat_type>;

//     using impl_type = impl::ImplicitStepperStandardImpl<
//       ind_var_type, state_type, residual_type, jacobian_type, policy_type>;
//     return impl::create_implicit_stepper_impl<
//       impl_type>(schemeName, policy_type(std::move(ws)));
//   }
//   else{
//     using policy_type = impl::ResidualJacobianStandardPolicy<
//       wrap_type, ind_var_type, state_type, residual_type, jacobian_type>;

//     using impl_type = impl::ImplicitStepperStandardImpl<
//       ind_var_type, state_type, residual_type, jacobian_type, policy_type>;
//     return impl::create_implicit_stepper_impl<
//       impl_type>(schemeName, policy_type(std::move(ws)));
//   }
// }

// template<class SystemType>
// auto create_implicit_stepper(StepScheme schemeName,
// 			     std::unique_ptr<SystemType> odeSystem)
// {
//   using system_type = mpl::remove_cvref_t<SystemType>;
//   static_assert(RealValuedOdeSystemFusingRhsAndJacobian<system_type>::value ||
// 		RealValuedCompleteOdeSystem<system_type>::value,
// 		"To create an implicit stepper you must provide an ODE system that meets either the RealValuedOdeSystemFusingRhsAndJacobian or the RealValuedCompleteOdeSystem concept.");
//   static constexpr bool complete_system = RealValuedCompleteOdeSystem<system_type>::value;

//   assert(
// 	 schemeName == StepScheme::BDF1 ||
// 	 schemeName == StepScheme::BDF2 ||
// 	 schemeName == StepScheme::CrankNicolson
// 	 );
//   if constexpr(complete_system){
//     // for systems with mass matrix CN not supported yet
//     assert(schemeName != StepScheme::CrankNicolson);
//   }

//   using ind_var_type  = typename system_type::independent_variable_type;
//   using state_type    = typename system_type::state_type;
//   using residual_type = typename system_type::rhs_type;
//   using jacobian_type = typename system_type::jacobian_type;

//   using wrap_type = impl::SystemInternalWrapper<1, system_type>;
//   wrap_type ws(std::move(odeSystem));

//   if constexpr (complete_system){
//     using mass_mat_type = typename system_type::mass_matrix_type;

//     using policy_type = impl::ResidualJacobianWithMassMatrixStandardPolicy<
//       wrap_type, ind_var_type, state_type,
//       residual_type, jacobian_type, mass_mat_type>;

//     using impl_type = impl::ImplicitStepperStandardImpl<
//       ind_var_type, state_type, residual_type, jacobian_type, policy_type>;
//     return impl::create_implicit_stepper_impl<
//       impl_type>(schemeName, policy_type(std::move(ws)));
//   }
//   else{
//     using policy_type = impl::ResidualJacobianStandardPolicy<
//       wrap_type, ind_var_type, state_type, residual_type, jacobian_type>;

//     using impl_type = impl::ImplicitStepperStandardImpl<
//       ind_var_type, state_type, residual_type, jacobian_type, policy_type>;
//     return impl::create_implicit_stepper_impl<
//       impl_type>(schemeName, policy_type(std::move(ws)));
//   }
// }

template<class ResidualJacobianPolicyType>
auto create_implicit_stepper_with_custom_policy(StepScheme schemeName,
						ResidualJacobianPolicyType && policy)
{
  using policy_type = mpl::remove_cvref_t<ResidualJacobianPolicyType>;
  static_assert(ImplicitResidualJacobianPolicy<policy_type>::value,
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


// /**
//  * Constructs and returns a stepper for doing implicit time integration on
//  * a user-defined system that conforms to Pressio's Fully Discrete
//  * API with Jacobian support.
//  *
//  * - TotalNumberOfDesiredStates:
//  *   This represents the number of discrete solution states used by the stepper.
//  *   Must be exactly 2 for now.
//  *
//  * - odeSystem:
//  *     A reference to the user-defined odeSystem.
//  *     Its type must conform to the RealValuedFullyDiscreteSystemWithJacobian concept.
//  *     This odeSystem is passed to the stepper and will be used during time integration
//  *     to evaluate the discrete residual and Jacobian.
//  *     Note that the odeSystem object must not go out of scope
//  *     while the stepper object is used.
//  *
//  *
//  * Static Assertions:
//  * - Ensures that TotalNumberOfDesiredStates is exactly 2.
//  * - odeSystem satisfies the RealValuedFullyDiscreteSystemWithJacobian concept
//  *
//  * Returns:
//  * - The stepper instance. Note that the type of this stepper is unspecified and
//  *   users should not rely on it, but should only rely on its interface.
//  *   This stepper object can be used to advance the odeSystem in
//  *   time using an implicit scheme.
//  *
//  */
// template<int TotalNumberOfDesiredStates, class SystemType>
// auto create_implicit_stepper(SystemType const & odeSystem)
// {
//   static_assert(TotalNumberOfDesiredStates == 2,
// 		"create_implicit_stepper currently only supports 2 total states");

//   using system_type = mpl::remove_cvref_t<SystemType>;
//   static_assert(RealValuedFullyDiscreteSystemWithJacobian<
// 		system_type, TotalNumberOfDesiredStates>::value,
// 		"The ode system passed does not meet the FullyDiscrete API");

//   using ind_var_type  = typename system_type::independent_variable_type;
//   using state_type    = typename system_type::state_type;
//   using residual_type = typename system_type::discrete_residual_type;
//   using jacobian_type = typename system_type::discrete_jacobian_type;

//   // in this case the user-provided ode system is captured by "reference"
//   // meaning that the user need to ensure the odeSystem lives long enough
//   // for the stepper to be used
//   using system_wrapper_type = impl::SystemInternalWrapper<0, system_type>;

//   using stepper_type = impl::StepperArbitrary<
//     TotalNumberOfDesiredStates, system_wrapper_type, ind_var_type,
//     state_type, residual_type, jacobian_type
//     >;
//   return stepper_type(system_wrapper_type(odeSystem));
// }

// template<int TotalNumberOfDesiredStates, class SystemType>
// auto create_implicit_stepper(std::unique_ptr<SystemType> odeSystem)
// {
//   static_assert(TotalNumberOfDesiredStates == 2,
// 		"create_implicit_stepper currently only supports 2 total states");

//   using system_type = mpl::remove_cvref_t<SystemType>;
//   static_assert(RealValuedFullyDiscreteSystemWithJacobian<
// 		system_type, TotalNumberOfDesiredStates>::value,
// 		"The ode system passed does not meet the FullyDiscrete API");

//   using ind_var_type  = typename system_type::independent_variable_type;
//   using state_type    = typename system_type::state_type;
//   using residual_type = typename system_type::discrete_residual_type;
//   using jacobian_type = typename system_type::discrete_jacobian_type;

//   using system_wrapper_type = impl::SystemInternalWrapper<1, system_type>;
//   system_wrapper_type sw(std::move(odeSystem));

//   using stepper_type = impl::StepperArbitrary<
//     TotalNumberOfDesiredStates, system_wrapper_type, ind_var_type,
//     state_type, residual_type, jacobian_type
//     >;
//   return stepper_type(std::move(sw));
// }

}} // end namespace pressio::ode
#endif  // PRESSIOROM_ODE_ODE_CREATE_IMPLICIT_STEPPER_HPP_
