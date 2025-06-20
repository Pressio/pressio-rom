/*
//@HEADER
// ************************************************************************
//
// ode_explicit_stepper.hpp
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

#ifndef PRESSIO_ODE_IMPL_ODE_EXPLICIT_STEPPER_WITHOUT_MASS_MATRIX_HPP_
#define PRESSIO_ODE_IMPL_ODE_EXPLICIT_STEPPER_WITHOUT_MASS_MATRIX_HPP_

#include <array>

namespace pressio{ namespace ode{ namespace impl{

// this class is NOT meant for direct instantiation.
// One needs to use the public create_* functions because
// templates are handled and passed properly there.
template<
  class SysWrapperType,
  class StateType,
  class IndVarType,
  class RightHandSideType
  >
class ExplicitStepperNoMassMatrixImpl{

public:
  using independent_variable_type  = IndVarType;
  using state_type  = StateType;

private:
  StepScheme name_;
  SysWrapperType systemObj_;
  std::vector<RightHandSideType> rhsInstances_;
  StateType auxiliaryState_;

public:
  ExplicitStepperNoMassMatrixImpl() = default;
  ExplicitStepperNoMassMatrixImpl(const ExplicitStepperNoMassMatrixImpl &) = delete;
  ExplicitStepperNoMassMatrixImpl & operator=(const ExplicitStepperNoMassMatrixImpl &) = delete;
  ExplicitStepperNoMassMatrixImpl(ExplicitStepperNoMassMatrixImpl &&) = default;
  ExplicitStepperNoMassMatrixImpl & operator=(ExplicitStepperNoMassMatrixImpl &&) = default;
  ~ExplicitStepperNoMassMatrixImpl() = default;

  ExplicitStepperNoMassMatrixImpl(ode::ForwardEuler,
				  SysWrapperType && sysObjW)
    : name_(StepScheme::ForwardEuler),
      systemObj_(std::move(sysObjW)),
      rhsInstances_{systemObj_.get().createRhs()},
      auxiliaryState_{systemObj_.get().createState()}
  {}

  ExplicitStepperNoMassMatrixImpl(ode::RungeKutta4,
				  SysWrapperType && sysObjW)
    : name_(StepScheme::RungeKutta4),
      systemObj_(std::move(sysObjW)),
      rhsInstances_{systemObj_.get().createRhs(),
		    systemObj_.get().createRhs(),
		    systemObj_.get().createRhs(),
		    systemObj_.get().createRhs()},
      auxiliaryState_{systemObj_.get().createState()}
  {}

  ExplicitStepperNoMassMatrixImpl(ode::AdamsBashforth2,
				  SysWrapperType && sysObjW)
    : name_(StepScheme::AdamsBashforth2),
      systemObj_(std::move(sysObjW)),
      rhsInstances_{systemObj_.get().createRhs(),
		    systemObj_.get().createRhs()},
      auxiliaryState_{systemObj_.get().createState()}
  {}

  ExplicitStepperNoMassMatrixImpl(ode::SSPRungeKutta3,
				  SysWrapperType && sysObjW)
    : name_(StepScheme::SSPRungeKutta3),
      systemObj_(std::move(sysObjW)),
      rhsInstances_{systemObj_.get().createRhs()},
      auxiliaryState_{systemObj_.get().createState()}
  {}

public:
  void operator()(StateType & odeState,
		  const ::pressio::ode::StepStartAt<independent_variable_type> & stepStartVal,
		  ::pressio::ode::StepCount step,
		  ::pressio::ode::StepSize<independent_variable_type> stepSize)
  {
    auto dummyRhsObserver = [](::pressio::ode::StepCount /*unused*/,
			       ::pressio::ode::IntermediateStepCount /*unused*/,
			       const independent_variable_type & /*unused*/,
			       const RightHandSideType & /*unused*/)
    { /*no op*/ };

    (*this)(odeState, stepStartVal, step, stepSize, dummyRhsObserver);
  }

  template<class RhsObserverType>
  void operator()(StateType & odeState,
		  const ::pressio::ode::StepStartAt<independent_variable_type> & stepStartVal,
		  ::pressio::ode::StepCount step,
		  ::pressio::ode::StepSize<independent_variable_type> stepSize,
		  RhsObserverType & rhsObserver)
  {
    if (name_ == ode::StepScheme::ForwardEuler){
      doStepImpl(ode::ForwardEuler(), odeState,
		 stepStartVal.get(), stepSize.get(),
		 step, rhsObserver);
    }

    else if (name_ == ode::StepScheme::RungeKutta4){
      doStepImpl(ode::RungeKutta4(), odeState,
		 stepStartVal.get(), stepSize.get(),
		 step, rhsObserver);
    }

    else if (name_ == ode::StepScheme::AdamsBashforth2){
      doStepImpl(ode::AdamsBashforth2(), odeState,
		 stepStartVal.get(), stepSize.get(),
		 step, rhsObserver);
    }

    else if (name_ == ode::StepScheme::SSPRungeKutta3){
      doStepImpl(ode::SSPRungeKutta3(), odeState,
		 stepStartVal.get(), stepSize.get(),
		 step, rhsObserver);
    }
  }

private:
  template<class RhsObserverType>
  void doStepImpl(ode::ForwardEuler,
		  StateType & odeState,
		  const independent_variable_type & stepStartTime,
		  const independent_variable_type & stepSize,
		  ::pressio::ode::StepCount stepNumber,
		  RhsObserverType & rhsObserver)
  {
    PRESSIOLOG_DEBUG("euler forward stepper: do step");

    //eval rhs
    auto & rhs = rhsInstances_[0];
    systemObj_.get().rhs(odeState, stepStartTime, rhs);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(0), stepStartTime, rhs);

    // y = y + stepSize * rhs
    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    constexpr auto one = static_cast<scalar_type>(1);
    ::pressio::ops::update(odeState, one, rhs, stepSize);
  }

  template<class RhsObserverType>
  void doStepImpl(ode::AdamsBashforth2,
		  StateType & odeState,
		  const independent_variable_type & stepStartTime,
		  const independent_variable_type & stepSize,
		  ::pressio::ode::StepCount stepNumber,
		  RhsObserverType & rhsObserver)
  {
    PRESSIOLOG_DEBUG("adams-bashforth2 stepper: do step");

    // // y_n+1 = y_n + stepSize*[ (3/2)*f(y_n, t_n) - (1/2)*f(y_n-1, t_n-1) ]

    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    using cnst = ::pressio::ode::constants::Constants<scalar_type>;
    const auto cfn   = cnst::threeOvTwo()*stepSize;
    const auto cfnm1 = cnst::negOneHalf()*stepSize;

    if (stepNumber.get()==1){
      // use Euler forward or we could use something else here maybe RK4
      auto & rhs = rhsInstances_[0];
      systemObj_.get().rhs(odeState, stepStartTime, rhs);
      rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(0), stepStartTime, rhs);
      ::pressio::ops::update(odeState, cnst::one(), rhs, stepSize);
    }
    else{
      auto & fn   = rhsInstances_[0];
      auto & fnm1 = rhsInstances_[1];
      // fn -> fnm1
      ::pressio::ops::deep_copy(fnm1, fn);

      systemObj_.get().rhs(odeState, stepStartTime, fn);
      rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(0), stepStartTime, fn);
      ::pressio::ops::update(odeState, cnst::one(), fn, cfn, fnm1, cfnm1);
    }
  }

  template<class RhsObserverType>
  void doStepImpl(ode::SSPRungeKutta3,
		  StateType & odeState,
		  const independent_variable_type & stepStartTime,
		  const independent_variable_type & stepSize,
		  ::pressio::ode::StepCount stepNumber,
		  RhsObserverType & rhsObserver)
  {
    PRESSIOLOG_DEBUG("ssprk3 stepper: do step");

    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    using cnst = ::pressio::ode::constants::Constants<scalar_type>;

    auto & rhs0 = rhsInstances_[0];

    // see e.g. https://gkeyll.readthedocs.io/en/latest/dev/ssp-rk.html

    const scalar_type stepSize_half{stepSize/cnst::two()};
    const independent_variable_type t_phalf{stepStartTime + stepSize_half};
    const independent_variable_type t_next{stepStartTime + stepSize};

    // rhs(u_n, t_n)
    systemObj_.get().rhs(odeState, stepStartTime, rhs0);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(0), stepStartTime, rhs0);
    // u_1 = u_n + stepSize * rhs(u_n, t_n)
    ::pressio::ops::update(auxiliaryState_, cnst::zero(),
                           odeState,        cnst::one(),
                           rhs0,            stepSize);

    // rhs(u_1, t_n+stepSize)
    systemObj_.get().rhs(auxiliaryState_, t_next, rhs0);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(1), t_next, rhs0);
    // u_2 = 3/4*u_n + 1/4*u_1 + 1/4*stepSize*rhs(u_1, t_n+stepSize)
    ::pressio::ops::update(
      auxiliaryState_, cnst::fourInv(), odeState,
      cnst::threeOvFour(), rhs0, cnst::fourInv()*stepSize
    );

    // rhs(u_2, t_n + 0.5*stepSize)
    systemObj_.get().rhs(auxiliaryState_, t_phalf, rhs0);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(2), t_phalf, rhs0);
    // u_n+1 = 1/3*u_n + 2/3*u_2 + 2/3*stepSize*rhs(u_2, t_n+0.5*stepSize)
    ::pressio::ops::update(
      odeState, cnst::oneOvThree(), auxiliaryState_,
      cnst::twoOvThree(), rhs0, cnst::twoOvThree()*stepSize
    );
  }

  template<class RhsObserverType>
  void doStepImpl(ode::RungeKutta4,
		  StateType & odeState,
		  const independent_variable_type & stepStartTime,
		  const independent_variable_type & stepSize,
		  ::pressio::ode::StepCount stepNumber,
		  RhsObserverType & rhsObserver)
  {
    PRESSIOLOG_DEBUG("rk4 stepper: do step");

    auto & rhs1 = rhsInstances_[0];
    auto & rhs2 = rhsInstances_[1];
    auto & rhs3 = rhsInstances_[2];
    auto & rhs4 = rhsInstances_[3];

    using scalar_type = typename ::pressio::Traits<StateType>::scalar_type;
    using cnst = ::pressio::ode::constants::Constants<scalar_type>;

    const scalar_type stepSize_half{stepSize/cnst::two()};
    const independent_variable_type t_phalf{stepStartTime + stepSize_half};
    const independent_variable_type t_next{stepStartTime + stepSize};
    const scalar_type stepSize6{stepSize/cnst::six()};
    const scalar_type stepSize3{stepSize/cnst::three()};

    // stage 1:
    // rhs1 = rhs(y_n, t_n)
    systemObj_.get().rhs(odeState, stepStartTime, rhs1);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(0), stepStartTime, rhs1);

    // stage 2:
    // ytmp = y + rhs1*stepSize_half;
    this->rk4_stage_update_impl(auxiliaryState_, odeState, rhs1, stepSize_half);
    // rhs2 = rhs(y_tmp, t_n+stepSize/2)
    systemObj_.get().rhs(auxiliaryState_, t_phalf, rhs2);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(1), t_phalf, rhs2);

    // stage 3:
    // ytmp = y + rhs2*stepSize_half;
    this->rk4_stage_update_impl(auxiliaryState_, odeState, rhs2, stepSize_half);
    // rhs3 = rhs(y_tmp)
    systemObj_.get().rhs(auxiliaryState_, t_phalf, rhs3);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(2), t_phalf, rhs3);

    // stage 4:
    // ytmp = y + rhs3*stepSize;
    this->rk4_stage_update_impl(auxiliaryState_, odeState, rhs3, stepSize);
    // rhs3 = rhs(y_tmp)
    systemObj_.get().rhs(auxiliaryState_, t_next, rhs4);
    rhsObserver(stepNumber, ::pressio::ode::IntermediateStepCount(3), t_next, rhs4);

    // y_n += stepSize/6 * ( rhs1 + 2*rhs2 + 2*rhs3 + rhs4 )
    ::pressio::ops::update(odeState, cnst::one(),
			   rhs1, stepSize6, rhs2, stepSize3,
			   rhs3, stepSize3, rhs4, stepSize6);
  }

  template<class rhs_t, class FactorType>
  void rk4_stage_update_impl(StateType & yIn,
			     const StateType & stateIn,
			     const rhs_t & rhsIn,
			     const FactorType & rhsFactor)
  {
    using scalar_type   = typename ::pressio::Traits<StateType>::scalar_type;
    using cnst = ::pressio::ode::constants::Constants<scalar_type>;
    ::pressio::ops::update(yIn, cnst::zero(), stateIn, cnst::one(), rhsIn, rhsFactor);
  }

};

}}}//end namespace pressio::ode::explicitmethods::impl
#endif  // PRESSIO_ODE_IMPL_ODE_EXPLICIT_STEPPER_WITHOUT_MASS_MATRIX_HPP_
