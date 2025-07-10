/*
//@HEADER
// ************************************************************************
//
// ode_implicit_discrete_residual.hpp
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

#ifndef PRESSIOROM_ODE_IMPL_ODE_IMPLICIT_DISCRETE_RESIDUAL_HPP_
#define PRESSIOROM_ODE_IMPL_ODE_IMPLICIT_DISCRETE_RESIDUAL_HPP_

namespace pressio{ namespace ode{ namespace impl{

/*
  BDF1 residual:
  R(y_n+1) = y_n+1 - y_n - dt*f(t_n+1, y_n+1)

  on entry R contains the application RHS: R = f(t_n+1, y_n+1, ...)
*/
template <
  class StateType,
  class ResidualType,
  class StencilStatesContainerType,
  class StepSizeType
  >
std::enable_if_t<
  ::pressio::all_have_traits_and_same_scalar<StateType, ResidualType>::value
  && std::is_convertible<
    StepSizeType, typename Traits<ResidualType>::scalar_type
    >::value
  >
discrete_residual(::pressio::ode::BDF1,
		  const StateType & y_np1,
		  ResidualType & R,
		  const StencilStatesContainerType & stencilStates,
		  const StepSizeType & dt)
{

  using sc_t = typename ::pressio::Traits<ResidualType>::scalar_type;
  const sc_t cf = ::pressio::ode::constants::bdf1<sc_t>::c_f_ * dt;
  constexpr sc_t cnp1 = ::pressio::ode::constants::bdf1<sc_t>::c_np1_;
  constexpr sc_t cn = ::pressio::ode::constants::bdf1<sc_t>::c_n_;
  const auto & y_n = stencilStates(::pressio::ode::n());

  ::pressio::ops::update(R, cf, y_np1, cnp1, y_n, cn);
}

/*
  BDF1 with MM residual:
  R(y_n+1) = M_n+1(y_n+1 - y_n) - dt*f(t_n+1, y_n+1)

  - on entry R is empty
  - on entry f contains the application f(t_n+1, y_n+1, ...)
  - on entry M contains mass matrix M_n+1
*/
template <
  class StateType,
  class MassMatrixType,
  class ResidualType,
  class StencilStatesContainerType,
  class StepSizeType
  >
std::enable_if_t<
  ::pressio::all_have_traits_and_same_scalar<StateType, ResidualType, MassMatrixType>::value
  && std::is_convertible<
    StepSizeType, typename Traits<ResidualType>::scalar_type
    >::value
  >
discrete_residual(::pressio::ode::BDF1,
		  const StateType & y_np1,
		  StateType & scratchState,
		  const ResidualType & f_np1,
		  const MassMatrixType & M_np1,
		  ResidualType & R,
		  const StencilStatesContainerType & stencilStates,
		  const StepSizeType & dt)
{

  using sc_t = typename ::pressio::Traits<ResidualType>::scalar_type;
  constexpr sc_t zero = static_cast<sc_t>(0);
  constexpr sc_t one  = static_cast<sc_t>(1);

  const sc_t cf	      = ::pressio::ode::constants::bdf1<sc_t>::c_f_ * dt;
  const auto & y_n = stencilStates(::pressio::ode::n());

  // scratchState = (y_n+1 - y_n)
  ::pressio::ops::update(scratchState, zero, y_np1, one, y_n, -one);
  // R = M_n+1 * scratchState
  ::pressio::ops::product(::pressio::nontranspose(), one, M_np1, scratchState, zero, R);
  // R = R - dt * f_n+1
  ::pressio::ops::update(R, one, f_np1, cf);
}

/*
  BDF2 residual: R(y_n+1) = y_n+1 - (4/3)*y_n + (1/3)*y_n-1 - (2/3)*dt*f(t_n+1, y_n+1)

  on entry R contains the application RHS: R = f(t_n+1, y_n+1, ...)
*/
template <
  class StateType,
  class ResidualType,
  class StencilStatesContainerType,
  class StepSizeType
  >
std::enable_if_t<
  ::pressio::all_have_traits_and_same_scalar<StateType, ResidualType>::value
  && std::is_convertible<
    StepSizeType, typename Traits<ResidualType>::scalar_type
    >::value
  >
discrete_residual(::pressio::ode::BDF2,
		  const StateType & y_np1,
		  ResidualType & R,
		  const StencilStatesContainerType & stencilStates,
		  const StepSizeType & dt)
{

  using sc_t = typename ::pressio::Traits<ResidualType>::scalar_type;

  constexpr sc_t cnp1 = ::pressio::ode::constants::bdf2<sc_t>::c_np1_;
  constexpr sc_t cn   = ::pressio::ode::constants::bdf2<sc_t>::c_n_;
  constexpr sc_t cnm1 = ::pressio::ode::constants::bdf2<sc_t>::c_nm1_;
  const sc_t cf	      = ::pressio::ode::constants::bdf2<sc_t>::c_f_ * dt;

  const auto & y_n = stencilStates(::pressio::ode::n());
  const auto & y_nm1 = stencilStates(::pressio::ode::nMinusOne());

  ::pressio::ops::update(R, cf, y_np1, cnp1, y_n, cn, y_nm1, cnm1);
}

/*
  BDF2 with MM:
   R(y_n+1) = M_n+1(y_n+1 - (4/3)*y_n + (1/3)*y_n-1) - (2/3)*dt*f(t_n+1, y_n+1)

  on entry R is empty
  on entry f_np1 contains f(t_n+1, y_n+1, ...)
  on output, R contains the discrete residual
*/
template <
  class StateType,
  class MassMatrixType,
  class ResidualType,
  class StencilStatesContainerType,
  class StepSizeType
  >
std::enable_if_t<
  ::pressio::all_have_traits_and_same_scalar<StateType, ResidualType, MassMatrixType>::value
  && std::is_convertible<
    StepSizeType, typename Traits<ResidualType>::scalar_type
    >::value
  >
discrete_residual(::pressio::ode::BDF2,
		  const StateType & y_np1,
		  StateType & scratchState,
		  const ResidualType & f_np1,
		  const MassMatrixType & M_np1,
		  ResidualType & R,
		  const StencilStatesContainerType & stencilStates,
		  const StepSizeType & dt)
{

  using sc_t = typename ::pressio::Traits<ResidualType>::scalar_type;
  using cnst = ::pressio::ode::constants::Constants<sc_t>;

  constexpr sc_t cnp1 = ::pressio::ode::constants::bdf2<sc_t>::c_np1_;
  constexpr sc_t cn   = ::pressio::ode::constants::bdf2<sc_t>::c_n_;
  constexpr sc_t cnm1 = ::pressio::ode::constants::bdf2<sc_t>::c_nm1_;
  const sc_t cf	      = ::pressio::ode::constants::bdf2<sc_t>::c_f_ * dt;

  const auto & y_n = stencilStates(::pressio::ode::n());
  const auto & y_nm1 = stencilStates(::pressio::ode::nMinusOne());

  // scratchState = y_n+1 - (4/3)*y_n + (1/3)*y_n-1
  ::pressio::ops::update(scratchState, cnst::zero(), y_np1, cnp1, y_n, cn, y_nm1, cnm1);
  // R = M_n+1 * scratchState
  ::pressio::ops::product(::pressio::nontranspose(), cnst::one(), M_np1, scratchState, cnst::zero(), R);
  // R = R - dt * f_n+1
  ::pressio::ops::update(R, cnst::one(), f_np1, cf);
}

/*
  CrankNicolson residual:
  R(y_n+1) = y_n+1 - y_n - 0.5*dt*[ f(t_n+1, y_n+1) + f(t_n, y_n) ]

  - on entry, R does not contain anything, should be fully overwritten
  - stencilStates contain: y_n
  - stencilVelocities contain f_n, f_n+1
*/
template <
  class StateType,
  class ResidualType,
  class StencilStatesContainerType,
  class StencilVelocitiesContainerType,
  class StepSizeType
  >
std::enable_if_t<
  ::pressio::all_have_traits_and_same_scalar<StateType, ResidualType>::value
  && std::is_convertible<
    StepSizeType, typename Traits<ResidualType>::scalar_type
    >::value
  >
discrete_residual(::pressio::ode::CrankNicolson,
		  const StateType & y_np1,
		  ResidualType & R,
		  const StencilStatesContainerType & stencilStates,
		  const StencilVelocitiesContainerType & stencilVelocities,
		  const StepSizeType & dt)
{

  using sc_t = typename ::pressio::Traits<ResidualType>::scalar_type;

  using cnst = ::pressio::ode::constants::cranknicolson<sc_t>;
  constexpr sc_t cnp1  = cnst::c_np1_;
  constexpr sc_t cn    = cnst::c_n_;
  constexpr sc_t cfn   = cnst::c_fn_;
  const sc_t cfnDt   = cfn*dt;
  const sc_t cfnp1Dt = cfn*dt;

  const auto & y_n   = stencilStates(::pressio::ode::n());
  const auto & f_n   = stencilVelocities(::pressio::ode::n());
  const auto & f_np1 = stencilVelocities(::pressio::ode::nPlusOne());

  ::pressio::ops::update(R, static_cast<sc_t>(0),
			 y_np1, cnp1,
			 y_n, cn,
			 f_n, cfnDt,
			 f_np1, cfnp1Dt);
}

}}}//end namespace pressio::ode::impl
#endif  // PRESSIOROM_ODE_IMPL_ODE_IMPLICIT_DISCRETE_RESIDUAL_HPP_
