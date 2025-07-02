/*
//@HEADER
// ************************************************************************
//
// ode_implicit_discrete_jacobian.hpp
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

#ifndef PRESSIOROM_ODE_IMPL_ODE_IMPLICIT_DISCRETE_JACOBIAN_HPP_
#define PRESSIOROM_ODE_IMPL_ODE_IMPLICIT_DISCRETE_JACOBIAN_HPP_

namespace pressio{ namespace ode{ namespace impl{

/*
  BDF1: J(y_n+1) = I - dt*df_n+1/dy_n+1
  on input jac contains  df_n+1/dy_n+1
  on output, jac contains the discrete jacobian
*/
template <class JacobianType, class StepSizeType>
std::enable_if_t<
  std::is_convertible<
    StepSizeType, typename Traits<JacobianType>::scalar_type
    >::value
  >
discrete_jacobian(::pressio::ode::BDF1,
		  JacobianType & jac,
		  const StepSizeType & dt)
{

  using sc_t = typename ::pressio::Traits<JacobianType>::scalar_type;
  constexpr sc_t cnp1   = ::pressio::ode::constants::bdf1<sc_t>::c_np1_;
  const sc_t cf   = ::pressio::ode::constants::bdf1<sc_t>::c_f_ * dt;
  ::pressio::ops::scale(jac, cf);
  ::pressio::ops::add_to_diagonal(jac, cnp1);
}

/*
  BDF1 WITH MM: J(y_n+1) = M_n+1 - dt*df_n+1/dy_n+1
  on input jac contains  df_n+1/dy_n+1
  on output, jac contains the discrete jacobian
*/
template <class JacobianType, class MassMatrixType, class StepSizeType>
std::enable_if_t<
  ::pressio::all_have_traits_and_same_scalar<JacobianType, MassMatrixType>::value
  && std::is_convertible<
    StepSizeType, typename Traits<JacobianType>::scalar_type
    >::value
  >
discrete_jacobian(::pressio::ode::BDF1,
		  JacobianType & jac,
		  const MassMatrixType & M_np1,
		  const StepSizeType & dt)
{

  using sc_t = typename ::pressio::Traits<JacobianType>::scalar_type;
  constexpr sc_t cnp1   = ::pressio::ode::constants::bdf1<sc_t>::c_np1_;
  const sc_t cf   = ::pressio::ode::constants::bdf1<sc_t>::c_f_ * dt;
  ::pressio::ops::update(jac, cf, M_np1, cnp1);
}

/*
  BDF2: J(y_n+1) = I - (2/3)*dt*df_n+1/dy_n+1
  - on input jac contains  df_n+1/dy_n+1
  - on output, jac contains the discrete jacobian
*/
template <class JacobianType, class StepSizeType>
std::enable_if_t<
  std::is_convertible<
    StepSizeType, typename Traits<JacobianType>::scalar_type
    >::value
  >
discrete_jacobian(::pressio::ode::BDF2,
		  JacobianType & jac,
		  const StepSizeType & dt)
{

  using sc_t = typename ::pressio::Traits<JacobianType>::scalar_type;
  constexpr sc_t cnp1   = ::pressio::ode::constants::bdf2<sc_t>::c_np1_;
  const sc_t cf   = ::pressio::ode::constants::bdf2<sc_t>::c_f_ * dt;
  ::pressio::ops::scale(jac, cf);
  ::pressio::ops::add_to_diagonal(jac, cnp1);
}

/*
  BDF2 WITH MM: J(y_n+1) = M_n+1 - (2/3)*dt*df_n+1/dy_n+1
  - on input jac contains  df_n+1/dy_n+1
  - on output, jac contains the discrete jacobian
*/
template <class JacobianType, class MassMatrixType, class StepSizeType>
std::enable_if_t<
  ::pressio::all_have_traits_and_same_scalar<JacobianType, MassMatrixType>::value
  && std::is_convertible<
    StepSizeType, typename Traits<JacobianType>::scalar_type
    >::value
  >
discrete_jacobian(::pressio::ode::BDF2,
		  JacobianType & jac,
		  const MassMatrixType & M_np1,
		  const StepSizeType & dt)
{

  using sc_t = typename ::pressio::Traits<JacobianType>::scalar_type;
  // constexpr sc_t one  = static_cast<scalar_type>(1);
  constexpr sc_t cnp1 = ::pressio::ode::constants::bdf2<sc_t>::c_np1_;
  const sc_t cf   = ::pressio::ode::constants::bdf2<sc_t>::c_f_ * dt;
  ::pressio::ops::update(jac, cf, M_np1, cnp1);
}

/*
  CRANK NICOLSON: J(y_n+1) = I - 0.5*dt*df_n+1/dy_n+1
  - on input jac contains  df_n+1/dy_n+1
  - on output, jac contains the discrete jacobian
*/
template <class JacobianType, class StepSizeType>
std::enable_if_t<
  std::is_convertible<
    StepSizeType, typename Traits<JacobianType>::scalar_type
    >::value
  >
discrete_jacobian(::pressio::ode::CrankNicolson,
		  JacobianType & jac,
		  const StepSizeType & dt)
{

  using sc_t = typename ::pressio::Traits<JacobianType>::scalar_type;
  using cnst = ::pressio::ode::constants::cranknicolson<sc_t>;
  constexpr sc_t cnp1  = cnst::c_np1_;
  const sc_t cf = cnst::c_fnp1_ * dt;
  ::pressio::ops::scale(jac, cf);
  ::pressio::ops::add_to_diagonal(jac, cnp1);
}

}}}//end namespace pressio::ode::impl
#endif  // PRESSIOROM_ODE_IMPL_ODE_IMPLICIT_DISCRETE_JACOBIAN_HPP_
