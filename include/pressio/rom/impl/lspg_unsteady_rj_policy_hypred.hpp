/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_policy_residual.hpp
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

#ifndef PRESSIOROM_ROM_IMPL_LSPG_UNSTEADY_RJ_POLICY_HYPRED_HPP_
#define PRESSIOROM_ROM_IMPL_LSPG_UNSTEADY_RJ_POLICY_HYPRED_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class IndVarType,
  class ReducedStateType,
  class LspgResidualType,
  class LspgJacobianType,
  class TrialSubspaceType,
  class FomSystemType,
  class HypRedUpdaterType
  >
class LspgUnsteadyResidualJacobianPolicyHypRed
{
  using fom_states_mgr_type = FomStatesManager<TrialSubspaceType>;

public:
  // required
  using independent_variable_type = IndVarType;
  using state_type    = ReducedStateType;
  using residual_type = LspgResidualType;
  using jacobian_type = LspgJacobianType;

public:
  LspgUnsteadyResidualJacobianPolicyHypRed(const TrialSubspaceType & trialSubspace,
					   const FomSystemType & fomSystem,
					   std::unique_ptr<fom_states_mgr_type> fomStatesManager,
					   const HypRedUpdaterType & hrUpdater)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      hypRedUpdater_(hrUpdater),
      fomStateHelperInstance_(trialSubspace.createFullState()),
      fomStatesManager_(std::move(fomStatesManager))
  {}

  LspgUnsteadyResidualJacobianPolicyHypRed(LspgUnsteadyResidualJacobianPolicyHypRed const &) = delete;
  LspgUnsteadyResidualJacobianPolicyHypRed& operator=(LspgUnsteadyResidualJacobianPolicyHypRed const&) = delete;
  LspgUnsteadyResidualJacobianPolicyHypRed(LspgUnsteadyResidualJacobianPolicyHypRed &&) = default;
  LspgUnsteadyResidualJacobianPolicyHypRed& operator=(LspgUnsteadyResidualJacobianPolicyHypRed &&) = default;

public:
  state_type createState() const{
    // this needs to create an instance of the reduced state
    return trialSubspace_.get().createReducedState();
  }

  residual_type createResidual() const{
    residual_type R(fomSystem_.get().createRhs());
    return R;
  }

  jacobian_type createJacobian() const{
    const auto phi = trialSubspace_.get().basisOfTranslatedSpace();
    return fomSystem_.get().createResultOfJacobianActionOn(phi);
  }

  template <class StencilStatesContainerType, class StencilRhsContainerType>
  void operator()(::pressio::ode::StepScheme odeSchemeName,
		  const state_type & predictedReducedState,
		  const StencilStatesContainerType & reducedStatesStencilManager,
		  StencilRhsContainerType & fomRhsStencilManger,
		  const ::pressio::ode::StepEndAt<IndVarType> & rhsEvaluationTime,
		  ::pressio::ode::StepCount step,
		  const ::pressio::ode::StepSize<IndVarType> & dt,
		  residual_type & R,
		  std::optional<jacobian_type *> Jo) const
  {

    if (odeSchemeName == ::pressio::ode::StepScheme::BDF1)
    {
      (*this).template compute_impl_bdf<ode::BDF1>
	(predictedReducedState, reducedStatesStencilManager,
	 rhsEvaluationTime.get(), dt.get(), step.get(), R, Jo);
    }

    else if (odeSchemeName == ::pressio::ode::StepScheme::BDF2)
    {
      if (step.get() == ::pressio::ode::first_step_value){
	(*this).template compute_impl_bdf<ode::BDF1>
	  (predictedReducedState, reducedStatesStencilManager,
	   rhsEvaluationTime.get(), dt.get(), step.get(), R, Jo);
      }
      else{
	(*this).template compute_impl_bdf<ode::BDF2>
	  (predictedReducedState, reducedStatesStencilManager,
	   rhsEvaluationTime.get(), dt.get(), step.get(), R, Jo);
      }
    }

    else{
      throw std::runtime_error("Invalid choice of StepScheme for hyp-red unstedy LSPG");
    }
  }

private:
  template <class OdeTag, class StencilStatesContainerType>
  void compute_impl_bdf(const state_type & predictedReducedState,
			 const StencilStatesContainerType & reducedStatesStencilManager,
			 const IndVarType & rhsEvaluationTime,
			 const IndVarType & dt,
			 const typename ::pressio::ode::StepCount::value_type & step,
			 residual_type & R,
			 std::optional<jacobian_type *> & Jo) const
  {
    static_assert( std::is_same<OdeTag, ode::BDF1>::value ||
		   std::is_same<OdeTag, ode::BDF2>::value, "");

    /* the FOM state for the prediction has to be always recomputed
       regardless of whether the currentStepNumber changes since
       we might be inside a subiteration of the non-linear solve
       where the time step does not change but the predicted state does
     */
    fomStatesManager_->reconstructAtWithoutStencilUpdate(predictedReducedState,
							      ::pressio::ode::nPlusOne());
    const auto & fomStateAt_np1 = (*fomStatesManager_)(::pressio::ode::nPlusOne());

    /* previous FOM states should only be recomputed when the time step changes.
       The method below does not recompute all previous states, but only
       recomputes the n-th state and updates/shifts back all the other
       FOM states stored. */
    if (stepTracker_ != step){
      const auto & lspgStateAt_n = reducedStatesStencilManager(::pressio::ode::n());
      fomStatesManager_->reconstructAtWithStencilUpdate(lspgStateAt_n,
							     ::pressio::ode::n());
      stepTracker_ = step;
    }

    if (std::is_same<OdeTag, ode::BDF1>::value){

      /* BDF1 residual : R(y_n+1) = cnp1*y_n+1 + cn*y_n + cf*f(t_n+1, y_n+1)
	 which we do follows:
	 1. R = f(t_n+1, y_n+1)
	 2. fomStateHelpInstance_ = cnp1*y_np1 + cn*y_n
	 3. call the hypRedUpdater to handle the rest
      */
      // step 1
      fomSystem_.get().rhs(fomStateAt_np1, rhsEvaluationTime, R);
      // step 2
      const auto & fomStateAt_n = (*fomStatesManager_)(::pressio::ode::n());
      using fom_state_type = typename FomSystemType::state_type;
      using sc_t = typename ::pressio::Traits<fom_state_type>::scalar_type;
      constexpr auto cnp1 = ::pressio::ode::constants::bdf1<sc_t>::c_np1_;
      constexpr auto cn   = ::pressio::ode::constants::bdf1<sc_t>::c_n_;
      ::pressio::ops::update(fomStateHelperInstance_, static_cast<sc_t>(0),
			     fomStateAt_np1, cnp1,
			     fomStateAt_n, cn);

      // step 3
      const auto cf = ::pressio::ode::constants::bdf1<sc_t>::c_f_ * dt;
      hypRedUpdater_.get().updateSampleMeshOperandWithStencilMeshOne
	(R, cf, fomStateHelperInstance_, static_cast<sc_t>(1));
    }
    else{

      /* BDF2 residual : R(y_n+1) = cnp1*y_n+1 + cn*y_n + cnm1*y_nm1 + cf*f(t_n+1, y_n+1)
	 which we do follows:
	 1. R = f(t_n+1, y_n+1)
	 2. fomStateHelpInstance_ = cnp1*y_np1 + cn*y_n + cnm1*y_nm1
	 3. call the hypRedUpdater to handle the rest
      */
      // step 1
      fomSystem_.get().rhs(fomStateAt_np1, rhsEvaluationTime, R);
      // step 2
      const auto & fomStateAt_n = (*fomStatesManager_)(::pressio::ode::n());
      const auto & fomStateAt_nm1 = (*fomStatesManager_)(::pressio::ode::nMinusOne());

      using fom_state_type = typename FomSystemType::state_type;
      using sc_t = typename ::pressio::Traits<fom_state_type>::scalar_type;
      constexpr auto cnp1 = ::pressio::ode::constants::bdf2<sc_t>::c_np1_;
      constexpr auto cn   = ::pressio::ode::constants::bdf2<sc_t>::c_n_;
      constexpr auto cnm1 = ::pressio::ode::constants::bdf2<sc_t>::c_nm1_;
      ::pressio::ops::update(fomStateHelperInstance_, static_cast<sc_t>(0),
			     fomStateAt_np1, cnp1,
			     fomStateAt_n, cn,
			     fomStateAt_nm1, cnm1);

      // step 3
      const auto cf = ::pressio::ode::constants::bdf2<sc_t>::c_f_ * dt;
      hypRedUpdater_.get().updateSampleMeshOperandWithStencilMeshOne
	(R, cf, fomStateHelperInstance_, static_cast<sc_t>(1));
    }

    // deal with jacobian if needed
    if (Jo){
      auto & J = *Jo.value();
      // lspgJac = decoderJac + dt*coeff*J*decoderJac
      // where J is the d(fomrhs)/dy and coeff depends on the scheme.

      // first, store J*phi into J
      const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
      fomSystem_.get().applyJacobian(fomStateAt_np1, phi, rhsEvaluationTime, J);

      using sc_t = typename ::pressio::Traits<
	typename TrialSubspaceType::basis_matrix_type>::scalar_type;
      IndVarType cf = {};
      if (std::is_same<OdeTag, ode::BDF1>::value){
	cf = dt * ::pressio::ode::constants::bdf1<sc_t>::c_f_;
      }
      else{
	cf = dt * ::pressio::ode::constants::bdf2<sc_t>::c_f_;
      }

      hypRedUpdater_.get().updateSampleMeshOperandWithStencilMeshOne(J, cf, phi, static_cast<sc_t>(1));
    }
  }

private:
  using raw_step_type = typename ::pressio::ode::StepCount::value_type;
  static_assert(std::is_signed<raw_step_type>::value, "");
  // storedStep is used to keep track of which step we are at.
  // used to decide if we need to update/recompute the previous FOM states or not.
  // To avoid recomputing previous FOM states if we are not in a new time step.
  mutable raw_step_type stepTracker_ = -1;

  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::reference_wrapper<const HypRedUpdaterType> hypRedUpdater_;
  mutable typename FomSystemType::state_type fomStateHelperInstance_;

protected:
  std::unique_ptr<fom_states_mgr_type> fomStatesManager_;
};

}}}

#endif  // PRESSIOROM_ROM_IMPL_LSPG_UNSTEADY_RJ_POLICY_HYPRED_HPP_
