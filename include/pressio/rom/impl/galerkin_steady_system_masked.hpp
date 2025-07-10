/*
//@HEADER
// ************************************************************************
//
// galerkin_steady_system_masked.hpp
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

#ifndef PRESSIOROM_ROM_IMPL_GALERKIN_STEADY_SYSTEM_MASKED_HPP_
#define PRESSIOROM_ROM_IMPL_GALERKIN_STEADY_SYSTEM_MASKED_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
masked Galerkin problem represents:

   hypredOp masker[fom_r(phi x)] = 0

- fom_r is the fom "residual"
- phi is the basis
- masker is the masker

From this we get a "reduced" residual/jacobian:
R = phi^T hypredOp masker(fom_r(phi x))
J = phi^T hypredOp masker(dfom_r/dx(phi x) phi)
*/
template <
  class ReducedStateType,
  class ReducedResidualType,
  class ReducedJacobianType,
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HypRedOpType
  >
class GalerkinSteadyMaskedSystem
{

  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;

  // deduce the unmasked types
  using unmasked_fom_residual_type = typename FomSystemType::residual_type;
  using unmasked_fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createResultOfJacobianActionOn
	     (std::declval<basis_matrix_type const &>()));

  // deduce the masked types
  using masked_fom_residual_type =
    decltype(std::declval<MaskerType const>().createResultOfMaskActionOn
	     (std::declval<unmasked_fom_residual_type const &>()));

  using masked_fom_jac_action_result_type =
    decltype(std::declval<MaskerType const>().createResultOfMaskActionOn
	     (std::declval<unmasked_fom_jac_action_result_type const &>()));

public:
  // required aliases
  using state_type    = ReducedStateType;
  using residual_type = ReducedResidualType;
  using jacobian_type = ReducedJacobianType;

  GalerkinSteadyMaskedSystem() = delete;
  GalerkinSteadyMaskedSystem(const TrialSubspaceType & trialSubspace,
			     const FomSystemType & fomSystem,
			     const MaskerType & masker,
			     const HypRedOpType & hyperReducer)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      hyperReducer_(hyperReducer),
      masker_(masker),
      unMaskedFomResidual_(fomSystem.createResidual()),
      unMaskedFomJacAction_(fomSystem.createResultOfJacobianActionOn(trialSubspace_.get().basisOfTranslatedSpace())),
      maskedFomResidual_(masker.createResultOfMaskActionOn(unMaskedFomResidual_)),
      maskedFomJacAction_(masker.createResultOfMaskActionOn(unMaskedFomJacAction_))
  {}

  GalerkinSteadyMaskedSystem(GalerkinSteadyMaskedSystem const &) = delete;
  GalerkinSteadyMaskedSystem& operator=(GalerkinSteadyMaskedSystem const&) = delete;
  GalerkinSteadyMaskedSystem(GalerkinSteadyMaskedSystem &&) = default;
  GalerkinSteadyMaskedSystem& operator=(GalerkinSteadyMaskedSystem &&) = default;

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  residual_type createResidual() const{
    return impl::CreateGalerkinRhs<residual_type>()(trialSubspace_.get().dimension());
  }

  jacobian_type createJacobian() const{
    return impl::CreateGalerkinJacobian<jacobian_type>()(trialSubspace_.get().dimension());
  }

  void residualAndJacobian(const state_type & reducedState,
			   residual_type & reducedResidual,
			   std::optional<jacobian_type*> reducedJacobian) const
  {
    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    if (reducedJacobian){
      auto ja = std::optional<unmasked_fom_jac_action_result_type*>(&unMaskedFomJacAction_);
      fomSystem_.get().residualAndJacobianAction(fomState_, unMaskedFomResidual_,
						 phi, ja);
    }
    else{
      fomSystem_.get().residualAndJacobianAction(fomState_, unMaskedFomResidual_,
						 phi, {});
    }

    // Apply masking to the residual: maskedResidual = masker(unmaskedResidual)
    masker_(unMaskedFomResidual_, maskedFomResidual_);
    // Apply hyper-reduction to get the reduced residual: reducedResidual = hypRed(maskedResidual)
    hyperReducer_(maskedFomResidual_, reducedResidual);

    // do similarly for the jacobian
    if (reducedJacobian){
      masker_(unMaskedFomJacAction_, maskedFomJacAction_);
      hyperReducer_(maskedFomJacAction_, *reducedJacobian.value());
    }
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HypRedOpType> hyperReducer_;
  std::reference_wrapper<const MaskerType> masker_;

  // UNMASKED fom R,J instances
  mutable unmasked_fom_residual_type unMaskedFomResidual_;
  mutable unmasked_fom_jac_action_result_type unMaskedFomJacAction_;

  // MASKED fom R,J instances
  mutable masked_fom_residual_type maskedFomResidual_;
  mutable masked_fom_jac_action_result_type maskedFomJacAction_;
};

}}}
#endif  // PRESSIOROM_ROM_IMPL_GALERKIN_STEADY_SYSTEM_MASKED_HPP_
