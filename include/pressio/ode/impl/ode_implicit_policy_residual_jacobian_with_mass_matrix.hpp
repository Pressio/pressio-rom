/*
//@HEADER
// ************************************************************************
//
// ode_implicit_policy_residual_jacobian_with_mass_matrix.hpp
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

#ifndef PRESSIOROM_ODE_IMPL_ODE_IMPLICIT_POLICY_RESIDUAL_JACOBIAN_WITH_MASS_MATRIX_HPP_
#define PRESSIOROM_ODE_IMPL_ODE_IMPLICIT_POLICY_RESIDUAL_JACOBIAN_WITH_MASS_MATRIX_HPP_

namespace pressio{ namespace ode{ namespace impl{

template<
  class SysWrapperType,
  class IndVarType,
  class StateType,
  class ResidualType,
  class JacobianType,
  class MassMatrixType
  >
class ResidualJacobianWithMassMatrixStandardPolicy
{
public:
  // required
  using independent_variable_type = IndVarType;
  using state_type    = StateType;
  using residual_type = ResidualType;
  using jacobian_type = JacobianType;

public:
  ResidualJacobianWithMassMatrixStandardPolicy() = default;

  explicit ResidualJacobianWithMassMatrixStandardPolicy(SysWrapperType && sysObjW)
    : systemObj_(std::move(sysObjW)),
      scratchState_(systemObj_.get().createState()),
      massMatrix_(systemObj_.get().createMassMatrix()),
      rhs_(systemObj_.get().createRhs())
  {}

  ResidualJacobianWithMassMatrixStandardPolicy(const ResidualJacobianWithMassMatrixStandardPolicy &) = delete;
  ResidualJacobianWithMassMatrixStandardPolicy & operator=(const ResidualJacobianWithMassMatrixStandardPolicy &) = delete;
  ResidualJacobianWithMassMatrixStandardPolicy(ResidualJacobianWithMassMatrixStandardPolicy &&) = default;
  ResidualJacobianWithMassMatrixStandardPolicy & operator=(ResidualJacobianWithMassMatrixStandardPolicy &&) = default;
  ~ResidualJacobianWithMassMatrixStandardPolicy() = default;

public:
  const auto & viewMassMatrix() const {
    return massMatrix_;
  }

  StateType createState() const{
    StateType result(systemObj_.get().createState());
    return result;
  }

  ResidualType createResidual() const{
    ResidualType R(systemObj_.get().createRhs());
    return R;
  }

  JacobianType createJacobian() const{
    JacobianType JJ(systemObj_.get().createJacobian());
    return JJ;
  }

  template <
    class StencilStatesContainerType,
    class StencilVelocitiesContainerType
    >
  void operator()(StepScheme name,
		  const StateType & predictedState,
		  const StencilStatesContainerType & stencilStatesManager,
		  StencilVelocitiesContainerType & stencilVelocities,
		  const ::pressio::ode::StepEndAt<IndVarType> & rhsEvaluationTime,
		  ::pressio::ode::StepCount step,
		  const ::pressio::ode::StepSize<IndVarType> & dt,
		  ResidualType & R,
		  std::optional<jacobian_type*> Jo) const
  {

    if (name == StepScheme::BDF1){
      (*this).template compute_impl_bdf1
	(predictedState, stencilStatesManager,
	 stencilVelocities, rhsEvaluationTime.get(),
	 dt.get(), step.get(), R, Jo);
    }

    else if (name == StepScheme::BDF2){
      (*this).template compute_impl_bdf2
	(predictedState, stencilStatesManager,
	 stencilVelocities, rhsEvaluationTime.get(),
	 dt.get(), step.get(), R, Jo);
    }

    else if (name == StepScheme::CrankNicolson){
      throw std::runtime_error("CrankNicolson with mass matrix not yet implemented");
    }
  }

private:
  template <
  class StencilStatesContainerType,
  class StencilVelocitiesContainerType,
  class StepType>
  void compute_impl_bdf1(const StateType & predictedState,
			 const StencilStatesContainerType & stencilStatesManager,
			 StencilVelocitiesContainerType & /*unused*/,
			 const IndVarType & evalTime,
			 const IndVarType & dt,
			 const StepType & step,
			 ResidualType & R,
			 std::optional<jacobian_type*> & Jo) const
  {

    try{
      stepTracker_ = step;
      systemObj_.get().massMatrixAndRhsAndJacobian(predictedState, evalTime,
						   massMatrix_, rhs_, Jo);
      discrete_residual(BDF1(), predictedState, scratchState_, rhs_,
			massMatrix_, R, stencilStatesManager, dt);

      if (Jo){
	auto & Jv = *(Jo.value());
	discrete_jacobian(BDF1(), Jv, massMatrix_, dt);
      }
    }
    catch (::pressio::eh::VelocityFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }
  }

  template <
  class StencilStatesContainerType,
  class StencilVelocitiesContainerType,
  class StepType>
  void compute_impl_bdf2(const StateType & predictedState,
			 const StencilStatesContainerType & stencilStatesManager,
			 StencilVelocitiesContainerType & /*unused*/,
			 const IndVarType & evalTime,
			 const IndVarType & dt,
			 const StepType & step,
			 ResidualType & R,
       std::optional<jacobian_type*> Jo) const
  {

    stepTracker_ = step;
    auto cond = [=](){ return (step == ::pressio::ode::first_step_value); };

    try{
      systemObj_.get().massMatrixAndRhsAndJacobian(predictedState, evalTime,
					      massMatrix_, rhs_, Jo);
      if (cond()){
	discrete_residual(BDF1(), predictedState, scratchState_, rhs_,
			  massMatrix_, R, stencilStatesManager, dt);

	if (Jo){
	  auto & Jv = *(Jo.value());
	  discrete_jacobian(BDF1(), Jv, massMatrix_, dt);
	}
      }
      else{
	discrete_residual(BDF2(), predictedState, scratchState_, rhs_,
			  massMatrix_, R, stencilStatesManager, dt);
	if (Jo){
	  auto & Jv = *(Jo.value());
	  discrete_jacobian(BDF2(), Jv, massMatrix_, dt);
	}
      }
    }
    catch (::pressio::eh::VelocityFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }
  }

private:
  SysWrapperType systemObj_;
  mutable int32_t stepTracker_ = -1;
  mutable StateType scratchState_;
  mutable MassMatrixType massMatrix_;
  mutable ResidualType rhs_;
};

}}}//end namespace pressio::ode::implicitmethods::policy
#endif  // PRESSIOROM_ODE_IMPL_ODE_IMPLICIT_POLICY_RESIDUAL_JACOBIAN_WITH_MASS_MATRIX_HPP_
