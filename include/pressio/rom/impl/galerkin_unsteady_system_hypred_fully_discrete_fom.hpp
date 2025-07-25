/*
//@HEADER
// ************************************************************************
//
// galerkin_unsteady_system_hypred_fully_discrete_fom.hpp
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

#ifndef PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_HYPRED_FULLY_DISCRETE_FOM_HPP_
#define PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_HYPRED_FULLY_DISCRETE_FOM_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  hat{R} := hyperReducer fom_R(t_n, y_n+1, t_n, ...)
  hat{J} := hyperReducer fom_J phi

- fom_R is the fom discrete time residual
- fom_J is the fom discrete time jacobian
- hyperReducer is the hyper-red operator
*/

template <
  std::size_t n,
  class IndVarType,
  class ReducedStateType,
  class ReducedResidualType,
  class ReducedJacobianType,
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType
  >
class GalerkinHypRedFullyDiscreteSystem
{
  using raw_step_type = typename ::pressio::ode::StepCount::value_type;
  static_assert(std::is_signed<raw_step_type>::value, "");

  // deduce from the fom object the type of result of
  // applying the Jacobian to the basis
  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createResultOfDiscreteTimeJacobianActionOn
	     (std::declval<basis_matrix_type const &>()) );

public:
  // required
  using independent_variable_type = IndVarType;
  using state_type	          = ReducedStateType;
  using discrete_residual_type    = ReducedResidualType;
  using discrete_jacobian_type    = ReducedJacobianType;

  GalerkinHypRedFullyDiscreteSystem(const GalerkinHypRedFullyDiscreteSystem &) = delete;
  GalerkinHypRedFullyDiscreteSystem & operator=(const GalerkinHypRedFullyDiscreteSystem &) = delete;
  GalerkinHypRedFullyDiscreteSystem(GalerkinHypRedFullyDiscreteSystem &&) = default;
  GalerkinHypRedFullyDiscreteSystem & operator=(GalerkinHypRedFullyDiscreteSystem &&) = default;
  ~GalerkinHypRedFullyDiscreteSystem() = default;

  GalerkinHypRedFullyDiscreteSystem(const TrialSubspaceType & trialSubspace,
				    const FomSystemType & fomSystem,
				    const HyperReducerType & hyperReducer)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      hyperReducer_(hyperReducer),
      fomStatesManager_(create_galerkin_fom_states_manager<n>(trialSubspace)),
      fomResidual_(fomSystem.createDiscreteTimeResidual()),
      fomJacAction_(fomSystem.createResultOfDiscreteTimeJacobianActionOn(trialSubspace_.get().basisOfTranslatedSpace()))
  {}

public:
  state_type createState() const{
    // this needs to instantiate the reduced state
    return trialSubspace_.get().createReducedState();
  }

  discrete_residual_type createDiscreteResidual() const{
    // this needs to instantiate the reduced residual
    return impl::CreateGalerkinRhs<discrete_residual_type>()(trialSubspace_.get().dimension());
  }

  discrete_jacobian_type createDiscreteJacobian() const{
    // this needs to instantiate the reduced jacobian
    return impl::CreateGalerkinJacobian<discrete_jacobian_type>()(trialSubspace_.get().dimension());
  }

  template<typename step_t, std::size_t _n = n>
  std::enable_if_t< (_n==2) >
  discreteResidualAndJacobian(const step_t & currentStepNumber,
                              const independent_variable_type & time_np1,
                              const independent_variable_type & dt,
                              discrete_residual_type & galerkinResidual,
                              std::optional<discrete_jacobian_type *> galerkinJacobian,
                              const state_type & galerkin_state_np1,
                              const state_type & galerkin_state_n) const
  {
    doFomStatesReconstruction(currentStepNumber, galerkin_state_np1, galerkin_state_n);
    const bool computeJacobian = bool(galerkinJacobian);

    const auto & ynp1 = fomStatesManager_(::pressio::ode::nPlusOne());
    const auto & yn   = fomStatesManager_(::pressio::ode::n());
    try
    {
      queryFomOperators(currentStepNumber, time_np1, dt, computeJacobian, ynp1, yn);
    }
    catch (::pressio::eh::DiscreteTimeResidualFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }

    computeReducedOperators(time_np1, galerkinResidual, galerkinJacobian, computeJacobian);
  }

  template<typename step_t, std::size_t _n = n>
  std::enable_if_t< (_n==3) >
  discreteResidualAndJacobian(const step_t & currentStepNumber,
                              const independent_variable_type & time_np1,
                              const independent_variable_type & dt,
                              discrete_residual_type & galerkinResidual,
                              std::optional<discrete_jacobian_type *> galerkinJacobian,
                              const state_type & galerkin_state_np1,
                              const state_type & galerkin_state_n,
                              const state_type & galerkin_state_nm1) const
  {
    doFomStatesReconstruction(currentStepNumber, galerkin_state_np1,
                              galerkin_state_n, galerkin_state_nm1);

    const auto & ynp1 = fomStatesManager_(::pressio::ode::nPlusOne());
    const auto & yn   = fomStatesManager_(::pressio::ode::n());
    const auto & ynm1 = fomStatesManager_(::pressio::ode::nMinusOne());
    const bool computeJacobian = bool(galerkinJacobian);

    try{
      queryFomOperators(currentStepNumber, time_np1, dt, computeJacobian, ynp1, yn, ynm1);
    }
    catch (::pressio::eh::DiscreteTimeResidualFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }

    computeReducedOperators(time_np1, galerkinResidual, galerkinJacobian, computeJacobian);
  }

private:
  template<typename step_t, class ...States>
  void queryFomOperators(const step_t & currentStepNumber,
                         const independent_variable_type & time_np1,
                         const independent_variable_type & dt,
                         bool computeJacobian,
                         States && ... states) const
  {
    const auto phi = trialSubspace_.get().basisOfTranslatedSpace();

    if (computeJacobian){
      using op_ja_t = std::optional<fom_jac_action_result_type*>;
      fomSystem_.get().discreteTimeResidualAndJacobianAction(currentStepNumber, time_np1,
                                                             dt, fomResidual_, phi,
                                                             op_ja_t{&fomJacAction_},
                                                             std::forward<States>(states)...);
    }
    else{
      fomSystem_.get().discreteTimeResidualAndJacobianAction(currentStepNumber, time_np1,
                                                             dt, fomResidual_, phi, {},
                                                             std::forward<States>(states)...);
    }
  }

  void computeReducedOperators(const independent_variable_type & time_np1,
                               discrete_residual_type & galerkinResidual,
                               std::optional<discrete_jacobian_type *> galerkinJacobian,
                               bool computeJacobian) const
  {
    hyperReducer_(fomResidual_, time_np1, galerkinResidual);
    if (computeJacobian){
      hyperReducer_(fomJacAction_, time_np1, *(galerkinJacobian.value()));
    }
  }

  void doFomStatesReconstruction(const int32_t & step_number,
                                 const state_type & galerkin_state_np1) const
  {
    fomStatesManager_.reconstructAtWithoutStencilUpdate(galerkin_state_np1,
                                                        ::pressio::ode::nPlusOne());
  }

  void doFomStatesReconstruction(const int32_t & step_number,
				 const state_type & galerkin_state_np1,
				 const state_type & galerkin_state_n) const
  {
    /* the FOM state corresponding to the new predicted state has to be
     * recomputed every time regardless of the time step chaning or not,
     *  since we might be inside a non-linear solve
     * where the time step does not change but this residual method
     * is called multiple times. */
    fomStatesManager_.reconstructAtWithoutStencilUpdate(galerkin_state_np1,
							::pressio::ode::nPlusOne());

    /* previous FOM states should only be recomputed when the time step changes.
     * The method below does not recompute all previous states, but only
     * recomputes the n-th state and updates/shifts back all the other
     * FOM states stored. */
    if (stepTracker_ != step_number){
      fomStatesManager_.reconstructAtWithStencilUpdate(galerkin_state_n,
						       ::pressio::ode::n());
      stepTracker_ = step_number;
    }
  }

  void doFomStatesReconstruction(const int32_t & step_number,
				 const state_type & galerkin_state_np1,
				 const state_type & galerkin_state_n,
				 const state_type & galerkin_state_nm1) const
  {
    (void)galerkin_state_nm1;
    doFomStatesReconstruction(step_number, galerkin_state_np1, galerkin_state_n);
  }

protected:
  // storedStep is used to keep track of which step we are at.
  // used to decide if we need to update/recompute the previous FOM states or not.
  // To avoid recomputing previous FOM states if we are not in a new time step.
  mutable raw_step_type stepTracker_ = -1;

  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::reference_wrapper<const HyperReducerType> hyperReducer_;
  mutable FomStatesManager<TrialSubspaceType> fomStatesManager_;
  mutable typename FomSystemType::discrete_residual_type fomResidual_;
  mutable fom_jac_action_result_type fomJacAction_;
};

}}} // end pressio::rom::impl
#endif  // PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_HYPRED_FULLY_DISCRETE_FOM_HPP_
