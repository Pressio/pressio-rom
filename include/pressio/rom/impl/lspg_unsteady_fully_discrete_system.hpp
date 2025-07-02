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

#ifndef PRESSIOROM_ROM_IMPL_LSPG_UNSTEADY_FULLY_DISCRETE_SYSTEM_HPP_
#define PRESSIOROM_ROM_IMPL_LSPG_UNSTEADY_FULLY_DISCRETE_SYSTEM_HPP_

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
  using fom_state_type = typename TrialSubspaceType::full_state_type;

  using fom_states_mgr_type = FomStatesManager<TrialSubspaceType>;

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
			  std::unique_ptr<fom_states_mgr_type> fomStatesManager)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomStatesManager_(std::move(fomStatesManager))
  {}

  LspgFullyDiscreteSystem(LspgFullyDiscreteSystem const &) = delete;
  LspgFullyDiscreteSystem& operator=(LspgFullyDiscreteSystem const&) = delete;
  LspgFullyDiscreteSystem(LspgFullyDiscreteSystem &&) = default;
  LspgFullyDiscreteSystem& operator=(LspgFullyDiscreteSystem &&) = default;

public:
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

  template<class StepIntType, std::size_t _n = n>
  std::enable_if_t< (_n == 2)>
  preStepHook(StepIntType stepNumber,
	      IndVarType time, IndVarType dt,
	      const state_type & lspg_state_np1,
	      const state_type & lspg_state_n) const
  {
    /* this method is called once before starting a step, so we need to update
       the proposed state and the previous state */
    fomStatesManager_->reconstructAtWithoutStencilUpdate(lspg_state_np1,
     							      ::pressio::ode::nPlusOne());
    fomStatesManager_->reconstructAtWithoutStencilUpdate(lspg_state_n,
							      ::pressio::ode::n());

    static constexpr bool hashook = ::pressio::ode::has_const_pre_step_hook_method<
      mpl::remove_cvref_t<FomSystemType>, _n, StepIntType, IndVarType, fom_state_type
      >::value;
    if constexpr(hashook){
      const auto & ynp1 = (*fomStatesManager_)(::pressio::ode::nPlusOne());
      const auto & yn   = (*fomStatesManager_)(::pressio::ode::n());
      fomSystem_.get().preStepHook(stepNumber, time, dt, ynp1, yn);
    }
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
    /* the FOM state corresponding to the new predicted state has to be
     * recomputed every time since this method is called at every update
     * of the solution inside a non-linear solve loop
     * Note that here we don't reconstruct at "n" because that state was
     * reconstructed inside the preStepHook and has not changed inside the
     * solver loop
     */
    fomStatesManager_->reconstructAtWithoutStencilUpdate(lspg_state_np1,
							      ::pressio::ode::nPlusOne());

    const auto & ynp1 = (*fomStatesManager_)(::pressio::ode::nPlusOne());
    const auto & yn   = (*fomStatesManager_)(::pressio::ode::n());
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

protected:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::unique_ptr<fom_states_mgr_type> fomStatesManager_;
};

}}}
#endif  // PRESSIOROM_ROM_IMPL_LSPG_UNSTEADY_FULLY_DISCRETE_SYSTEM_HPP_
