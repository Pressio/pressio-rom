/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_discrete_time_default_system.hpp
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

#ifndef PRESSIO_ROM_IMPL_LSPG_UNSTEADY_FULLY_DISCRETE_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_LSPG_UNSTEADY_FULLY_DISCRETE_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  std::size_t n,
  class IndVarType,
  class ReducedStateType,
  class LspgResidualType,
  class LspgJacobianType,
  class TrialSubspaceType,
  class FomSystemType
  >
class LspgFullyDiscreteSystem
{
  using raw_step_type = typename ::pressio::ode::StepCount::value_type;
  static_assert(std::is_signed<raw_step_type>::value, "");

public:
  // required
  using independent_variable_type   = IndVarType;
  using state_type	            = ReducedStateType;
  using discrete_residual_type = LspgResidualType;
  using discrete_jacobian_type = LspgJacobianType;

  // needed by impl
  using residual_type = discrete_residual_type;
  using jacobian_type = discrete_jacobian_type;

  LspgFullyDiscreteSystem(const TrialSubspaceType & trialSubspace,
			  const FomSystemType & fomSystem,
			  LspgFomStatesManager<TrialSubspaceType> & fomStatesManager)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomStatesManager_(fomStatesManager)
  {}

  state_type createState() const{
    // this needs to create an instance of the reduced state
    return trialSubspace_.get().createReducedState();
  }

  discrete_residual_type createDiscreteResidual() const{
    discrete_residual_type R(fomSystem_.get().createDiscreteTimeResidual());
    return R;
  }

  discrete_jacobian_type createDiscreteJacobian() const{
    const auto phi = trialSubspace_.get().basisOfTranslatedSpace();
    discrete_jacobian_type J(fomSystem_.get().createResultOfDiscreteTimeJacobianActionOn(phi));
    return J;
  }

  template<class StepIntType>
  std::enable_if_t<
    ::pressio::ode::has_const_pre_step_hook_method<
      mpl::remove_cvref_t<FomSystemType>, StepIntType, IndVarType
      >::value
    >
  preStepHook(StepIntType stepNumber, IndVarType time, IndVarType dt) const{
    fomSystem_.get().preStepHook(stepNumber, time, dt);
  }

  template<typename step_t, std::size_t _n = n>
  std::enable_if_t< (_n==2) >
  discreteResidualAndJacobian(const step_t & currentStepNumber,
			      const independent_variable_type & time_np1,
			      const independent_variable_type & dt,
			      discrete_residual_type & R,
			      std::optional<discrete_jacobian_type*> Jo,
			      const state_type & lspg_state_np1,
			      const state_type & lspg_state_n) const
  {
    doFomStatesReconstruction(currentStepNumber, lspg_state_np1, lspg_state_n);
    const auto & ynp1 = fomStatesManager_(::pressio::ode::nPlusOne());
    const auto & yn   = fomStatesManager_(::pressio::ode::n());
    const auto phi = trialSubspace_.get().basisOfTranslatedSpace();

    try
    {
      fomSystem_.get().discreteTimeResidualAndJacobianAction(currentStepNumber, time_np1, dt,
							     R, phi, Jo, ynp1, yn);
    }
    catch (::pressio::eh::DiscreteTimeResidualFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }
  }

  template<typename step_t, std::size_t _n = n>
  std::enable_if_t< (_n==3) >
  discreteResidualAndJacobian(const step_t & currentStepNumber,
			      const independent_variable_type & time_np1,
			      const independent_variable_type & dt,
			      discrete_residual_type & R,
			      std::optional<discrete_jacobian_type*> Jo,
			      const state_type & lspg_state_np1,
			      const state_type & lspg_state_n,
			      const state_type & lspg_state_nm1) const
  {
    doFomStatesReconstruction(currentStepNumber, lspg_state_np1,
			      lspg_state_n, lspg_state_nm1);
    const auto & ynp1 = fomStatesManager_(::pressio::ode::nPlusOne());
    const auto & yn   = fomStatesManager_(::pressio::ode::n());
    const auto & ynm1 = fomStatesManager_(::pressio::ode::nMinusOne());
    const auto phi = trialSubspace_.get().basisOfTranslatedSpace();

    try{
      fomSystem_.get().discreteTimeResidualAndJacobianAction(currentStepNumber, time_np1, dt,
							     R, phi, Jo, ynp1, yn, ynm1);
    }
    catch (::pressio::eh::DiscreteTimeResidualFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }
  }

private:
  void doFomStatesReconstruction(const int32_t & step_number,
				 const state_type & lspg_state_np1) const
  {
    fomStatesManager_.get().reconstructAtWithoutStencilUpdate(lspg_state_np1,
							   ::pressio::ode::nPlusOne());
  }

  void doFomStatesReconstruction(const int32_t & step_number,
				 const state_type & lspg_state_np1,
				 const state_type & lspg_state_n) const
  {
    /* the FOM state corresponding to the new predicted state has to be
     * recomputed every time regardless of the time step chaning or not,
     *  since we might be inside a non-linear solve
     * where the time step does not change but this residual method
     * is called multiple times. */
    fomStatesManager_.get().reconstructAtWithoutStencilUpdate(lspg_state_np1,
							   ::pressio::ode::nPlusOne());

    /* previous FOM states should only be recomputed when the time step changes.
     * The method below does not recompute all previous states, but only
     * recomputes the n-th state and updates/shifts back all the other
     * FOM states stored. */
    if (stepTracker_ != step_number){
      fomStatesManager_.get().reconstructAtWithStencilUpdate(lspg_state_n,
							  ::pressio::ode::n());
      stepTracker_ = step_number;
    }
  }

  void doFomStatesReconstruction(const int32_t & step_number,
				 const state_type & lspg_state_np1,
				 const state_type & lspg_state_n,
				 const state_type & lspg_state_nm1) const
  {
    (void)lspg_state_nm1;
    doFomStatesReconstruction(step_number, lspg_state_np1, lspg_state_n);
  }

protected:
  // storedStep is used to keep track of which step we are at.
  // used to decide if we need to update/recompute the previous FOM states or not.
  // To avoid recomputing previous FOM states if we are not in a new time step.
  mutable raw_step_type stepTracker_ = -1;

  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::reference_wrapper<LspgFomStatesManager<TrialSubspaceType>> fomStatesManager_;
};

}}}
#endif  // PRESSIO_ROM_IMPL_LSPG_UNSTEADY_FULLY_DISCRETE_SYSTEM_HPP_
