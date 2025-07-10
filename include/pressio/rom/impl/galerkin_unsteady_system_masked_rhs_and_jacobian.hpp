/*
//@HEADER
// ************************************************************************
//
// galerkin_unsteady_system_masked_rhs_and_jacobian.hpp
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

#ifndef PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_AND_JACOBIAN_HPP_
#define PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_AND_JACOBIAN_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  masked implicit galerkin system represents:

     d hat{y}/dt = hyperReducer masker ( fom_rhs(phi*hat{y}, ...) )

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- hyperReducer is the hypred operator
so that it boils down to:

rhs = hyperReducer masked(fom_rhs(phi*hat{y}, ...))
rhs_jacobian = hyperReducer masked(d(fom_rhs(phi*hat{y}, ...))/dy phi)

*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class ReducedJacobianType,
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType
  >
class GalerkinMaskedOdeSystemRhsAndJacobian
{
  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;

  // deduce the unmasked types
  using unmasked_fom_rhs_type = typename FomSystemType::rhs_type;
  using unmasked_fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createResultOfJacobianActionOn
	     (std::declval<basis_matrix_type const &>()));

  // deduce the masked types
  // deduce the masked types
  using masked_fom_rhs_type =
    decltype(std::declval<MaskerType const>().createResultOfMaskActionOn
	     (std::declval<unmasked_fom_rhs_type const &>()));

  using masked_fom_jac_action_result_type =
    decltype(std::declval<MaskerType const>().createResultOfMaskActionOn
	     (std::declval<unmasked_fom_jac_action_result_type const &>()));

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using rhs_type      = ReducedRhsType;
  using jacobian_type = ReducedJacobianType;

  GalerkinMaskedOdeSystemRhsAndJacobian(const TrialSubspaceType & trialSubspace,
					const FomSystemType & fomSystem,
					const MaskerType & masker,
					const HyperReducerType & hyperReducer)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      hyperReducer_(hyperReducer),
      masker_(masker),
      unMaskedFomRhs_(fomSystem.createRhs()),
      unMaskedFomJacAction_(fomSystem.createResultOfJacobianActionOn(trialSubspace_.get().basisOfTranslatedSpace())),
      maskedFomRhs_(masker.createResultOfMaskActionOn(unMaskedFomRhs_)),
      maskedFomJacAction_(masker.createResultOfMaskActionOn(unMaskedFomJacAction_))
  {}

  GalerkinMaskedOdeSystemRhsAndJacobian(GalerkinMaskedOdeSystemRhsAndJacobian const &) = delete;
  GalerkinMaskedOdeSystemRhsAndJacobian& operator=(GalerkinMaskedOdeSystemRhsAndJacobian const&) = delete;
  GalerkinMaskedOdeSystemRhsAndJacobian(GalerkinMaskedOdeSystemRhsAndJacobian &&) = default;
  GalerkinMaskedOdeSystemRhsAndJacobian& operator=(GalerkinMaskedOdeSystemRhsAndJacobian &&) = default;

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  rhs_type createRhs() const{
    return impl::CreateGalerkinRhs<rhs_type>()(trialSubspace_.get().dimension());
  }

  jacobian_type createJacobian() const{
    return impl::CreateGalerkinJacobian<jacobian_type>()(trialSubspace_.get().dimension());
  }

  void rhsAndJacobian(const state_type & reducedState,
		      const IndVarType & rhsEvaluationTime,
		      rhs_type & reducedRhs,
		      std::optional<jacobian_type*> reducedJacobian) const
  {

    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);
    fomSystem_.get().rhs(fomState_, rhsEvaluationTime, unMaskedFomRhs_);
    masker_(unMaskedFomRhs_, maskedFomRhs_);
    hyperReducer_(maskedFomRhs_, rhsEvaluationTime, reducedRhs);

    if (reducedJacobian){
      const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
      fomSystem_.get().applyJacobian(fomState_, phi, rhsEvaluationTime, unMaskedFomJacAction_);
      masker_(unMaskedFomJacAction_, maskedFomJacAction_);
      hyperReducer_(maskedFomJacAction_, rhsEvaluationTime, *reducedJacobian.value());
    }
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReducerType> hyperReducer_;
  std::reference_wrapper<const MaskerType> masker_;

  // UNMASKED objects
  mutable unmasked_fom_rhs_type unMaskedFomRhs_;
  mutable unmasked_fom_jac_action_result_type unMaskedFomJacAction_;

  // MASKED objects
  mutable masked_fom_rhs_type maskedFomRhs_;
  mutable masked_fom_jac_action_result_type maskedFomJacAction_;
};

}}} // end pressio::rom::impl
#endif  // PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_AND_JACOBIAN_HPP_
