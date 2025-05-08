/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_arbitrary.hpp
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

#ifndef PRESSIO_ODE_IMPL_ODE_IMPLICIT_STEPPER_ARBITRARY_HPP_
#define PRESSIO_ODE_IMPL_ODE_IMPLICIT_STEPPER_ARBITRARY_HPP_

namespace pressio{ namespace ode{ namespace impl{

template<
  int n_states,
  class SystemType,
  class IndVarType,
  class StateType,
  class ResidualType,
  class JacobianType
  >
class StepperArbitrary
{
  static_assert(n_states == 2,
		"StepperArbitrary currently only supports 2 total states");

public:
  // required
  using independent_variable_type = IndVarType;
  using state_type	= StateType;
  using residual_type	= ResidualType;
  using jacobian_type	= JacobianType;

  // numAuxStates is the number of **auxiliary** states needed,
  // so total states - 1
  static constexpr std::size_t numAuxStates = n_states - 1;
  using tag_name = ::pressio::ode::ImplicitArbitrary;
  using stencil_states_t = ImplicitStencilStatesStaticContainer<StateType, numAuxStates>;

private:
  typename StepCount::value_type stepNumber_  = {};
  IndVarType dt_ = {};
  // for implicit integration, the time at the start of a step is
  // not equal to the time needed to evalaute the residual/jacobian.
  // so we keep track of both
  IndVarType timeAtStepStart_  = {};
  IndVarType rhsEvaluationTime_  = {};

  ::pressio::nonlinearsolvers::impl::InstanceOrReferenceWrapper<SystemType> systemObj_;
  stencil_states_t stencilStates_;
  // state object to ensure the strong guarantee for handling excpetions
  StateType recoveryState_;

public:
  StepperArbitrary() = delete;
  StepperArbitrary(const StepperArbitrary & other)  = default;
  StepperArbitrary & operator=(const StepperArbitrary & other) = delete;
  ~StepperArbitrary() = default;

  StepperArbitrary(SystemType && systemObj)
    : systemObj_(std::forward<SystemType>(systemObj)),
      stencilStates_(systemObj.createState()), //stencilstates handles right semantics
      recoveryState_{systemObj.createState()}
    {}

public:
  auto createState() const{
    return systemObj_.get().createState();
  }

  residual_type createResidual() const{
    return systemObj_.get().createDiscreteResidual();
  }

  jacobian_type createJacobian() const{
    return systemObj_.get().createDiscreteJacobian();
  }

  template<class SolverType, class ...Args>
  void operator()(StateType & odeState,
		  const ::pressio::ode::StepStartAt<independent_variable_type> & stepStartVal,
		  ::pressio::ode::StepCount stepNumber,
		  const ::pressio::ode::StepSize<independent_variable_type> & stepSize,
		  SolverType & solver,
		  Args && ...args)
  {
    PRESSIOLOG_DEBUG("arbitrary stepper: do step");

    // cache things before starting the step
    stepNumber_ = stepNumber.get();
    dt_ = stepSize.get();
    timeAtStepStart_ = stepStartVal.get();
    rhsEvaluationTime_ = stepStartVal.get() + dt_;

    // need to update auxiliary things **before** starting the step
    updateAuxiliaryStorage<numAuxStates>(odeState);

    callPreStepHookIfApplicable(odeState);

    try{
      solver.solve(*this, odeState, std::forward<Args>(args)...);
    }
    catch (::pressio::eh::NonlinearSolveFailure const & e)
    {
      rollBackStates<numAuxStates>(odeState);
      throw ::pressio::eh::TimeStepFailure();
    }
  }

  // 1 aux states, 2 total states
  template< std::size_t _numAuxStates = numAuxStates>
  std::enable_if_t< _numAuxStates==1 >
  residualAndJacobian(const state_type & odeState,
		      residual_type & R,
		      std::optional<jacobian_type*> Jo) const
  {
    const auto & yn = stencilStates_(ode::n());

    try{
      systemObj_.get().discreteResidualAndJacobian
	(stepNumber_, rhsEvaluationTime_, dt_, R, Jo, odeState, yn);
    }
    catch (::pressio::eh::DiscreteTimeResidualFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }
  }

private:
  void callPreStepHookIfApplicable(const StateType & odeState)
  {
    if constexpr (numAuxStates == 1 && has_const_pre_step_hook_method<
		  mpl::remove_cvref_t<SystemType>, n_states,
		  typename StepCount::value_type, IndVarType, state_type
		  >::value)
    {
      const auto & yn = stencilStates_(ode::n());
      systemObj_.get().preStepHook(stepNumber_, timeAtStepStart_, dt_,
				   odeState, yn);
    }
  }

  template<std::size_t nAux>
  std::enable_if_t<nAux==1>
  updateAuxiliaryStorage(const StateType & odeState)
  {
    auto & y_n = stencilStates_(ode::n());
    ::pressio::ops::deep_copy(recoveryState_, y_n);
    ::pressio::ops::deep_copy(y_n, odeState);
  }

  template<std::size_t nAux>
  std::enable_if_t<nAux==1>
  rollBackStates(StateType & odeState)
  {
    auto & y_n = stencilStates_(ode::n());
    ::pressio::ops::deep_copy(odeState, y_n);
    ::pressio::ops::deep_copy(y_n, recoveryState_);
  }
};

}}}
#endif  // PRESSIO_ODE_IMPL_ODE_IMPLICIT_STEPPER_ARBITRARY_HPP_
