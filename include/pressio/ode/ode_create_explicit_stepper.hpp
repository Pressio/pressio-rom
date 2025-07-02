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

#ifndef PRESSIO_ODE_ODE_CREATE_EXPLICIT_STEPPER_HPP_
#define PRESSIO_ODE_ODE_CREATE_EXPLICIT_STEPPER_HPP_

#include "./impl/ode_user_system_wrapper.hpp"
#include "./impl/ode_explicit_stepper_without_mass_matrix.hpp"
#include "./impl/ode_explicit_stepper_with_mass_matrix.hpp"
#include "./impl/ode_explicit_create_impl.hpp"

namespace pressio{ namespace ode{

template<class SystemType>
auto create_explicit_stepper(StepScheme schemeName,
			     SystemType const & odeSystem)
{
  using sys_type = mpl::remove_cvref_t<SystemType>;
  static_assert(
    RealValuedOdeSystem<sys_type>::value ||
    RealValuedOdeSystemFusingMassMatrixAndRhs<sys_type>::value,
    "To create an explicit stepper you must provide an ODE system that meets either the RealValuedOdeSystem or the RealValuedOdeSystemFusingMassMatrixAndRhs concept.");

  using ind_var_type = typename sys_type::independent_variable_type;
  using state_type   = typename sys_type::state_type;
  using rhs_type = typename sys_type::rhs_type;

  // in this case the user-provided ode system is captured by "reference"
  // meaning that the user need to ensure the odeSystem lives long enough
  // for the stepper to be used
  using system_wrapper_type = impl::SystemInternalWrapper<0, sys_type>;
  system_wrapper_type sw(odeSystem);

  static constexpr bool need_mass_matrix =
    RealValuedOdeSystemFusingMassMatrixAndRhs<sys_type>::value;
  if constexpr(need_mass_matrix){
    using mass_matrix_type = typename sys_type::mass_matrix_type;

    using impl_type = impl::ExplicitStepperWithMassMatrixImpl<
      system_wrapper_type, state_type, ind_var_type, rhs_type, mass_matrix_type>;
    return impl::create_explicit_stepper<impl_type>(schemeName, std::move(sw));
  }
  else{
    using impl_type = impl::ExplicitStepperNoMassMatrixImpl<
      system_wrapper_type, state_type, ind_var_type, rhs_type>;
    return impl::create_explicit_stepper<impl_type>(schemeName, std::move(sw));
  }
}

template<class SystemType>
auto create_explicit_stepper(StepScheme schemeName,
			     std::unique_ptr<SystemType> odeSystem)
{
  using sys_type = mpl::remove_cvref_t<SystemType>;
  static_assert(
    RealValuedOdeSystem<sys_type>::value ||
    RealValuedOdeSystemFusingMassMatrixAndRhs<sys_type>::value,
    "To create an explicit stepper you must provide an ODE system that meets either the RealValuedOdeSystem or the RealValuedOdeSystemFusingMassMatrixAndRhs concept.");

  using ind_var_type = typename sys_type::independent_variable_type;
  using state_type   = typename sys_type::state_type;
  using rhs_type = typename sys_type::rhs_type;

  // the user provides its ode system as a unique_ptr, so here we
  // need to take ownership of it and store it properly inside.
  using system_wrapper_type = impl::SystemInternalWrapper<1, sys_type>;
  system_wrapper_type sw(std::move(odeSystem));

  static constexpr bool need_mass_matrix =
    RealValuedOdeSystemFusingMassMatrixAndRhs<sys_type>::value;
  if constexpr(need_mass_matrix){
    using mass_matrix_type = typename sys_type::mass_matrix_type;

    using impl_type = impl::ExplicitStepperWithMassMatrixImpl<
      system_wrapper_type, state_type, ind_var_type, rhs_type, mass_matrix_type>;
    return impl::create_explicit_stepper<impl_type>(schemeName, std::move(sw));
  }
  else{
    using impl_type = impl::ExplicitStepperNoMassMatrixImpl<
      system_wrapper_type, state_type, ind_var_type, rhs_type>;
    return impl::create_explicit_stepper<impl_type>(schemeName, std::move(sw));
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
#endif  // PRESSIO_ODE_ODE_CREATE_EXPLICIT_STEPPER_HPP_
