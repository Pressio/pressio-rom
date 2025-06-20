
#ifndef PRESSIO_ROM_IMPL_LSPG_STEADY_SYSTEM_MASKED_HPP_
#define PRESSIO_ROM_IMPL_LSPG_STEADY_SYSTEM_MASKED_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
LSPG steady masked represents:

  min_x || mask[fom_r(phi x)]||

- fom_r is the fom "residual"
- phi is the basis
*/
template <
  class ReducedStateType,
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class PossiblyRefWrapperOperatorScalerType
  >
class LspgSteadyMaskedSystem
{

  using unmasked_fom_residual_type = typename FomSystemType::residual_type;
  using unmasked_fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createResultOfJacobianActionOn
      (std::declval<typename TrialSubspaceType::basis_matrix_type const &>())
      );

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
  using residual_type = masked_fom_residual_type;
  using jacobian_type = masked_fom_jac_action_result_type;

  // here _RawScalerType must be a template because it is the raw scaler type
  // which can be forwarded to the reference wrapper potentially
  template<class _RawScalerType>
  LspgSteadyMaskedSystem(const TrialSubspaceType & trialSubspace,
			 const FomSystemType & fomSystem,
			 const MaskerType & masker,
			 _RawScalerType && scaler)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      masker_(masker),
      unMaskedFomResidual_(fomSystem.createResidual()),
      unMaskedFomJacAction_(fomSystem.createResultOfJacobianActionOn(trialSubspace_.get().basisOfTranslatedSpace())),
      scaler_(std::forward<_RawScalerType>(scaler))
  {}

  LspgSteadyMaskedSystem(LspgSteadyMaskedSystem const &) = delete;
  LspgSteadyMaskedSystem& operator=(LspgSteadyMaskedSystem const&) = delete;
  LspgSteadyMaskedSystem(LspgSteadyMaskedSystem &&) = default;
  LspgSteadyMaskedSystem& operator=(LspgSteadyMaskedSystem &&) = default;

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  residual_type createResidual() const{
    auto tmp = fomSystem_.get().createResidual();
    return masker_.get().createResultOfMaskActionOn(tmp);
  }

  jacobian_type createJacobian() const{
    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    auto tmp = fomSystem_.get().createResultOfJacobianActionOn(phi);
    return masker_.get().createResultOfMaskActionOn(tmp);
  }

  void residualAndJacobian(const state_type & reducedState,
			   residual_type & lspgResidual,
			   std::optional<jacobian_type *> lspgJacobian) const
  {
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();

    if (lspgJacobian){
      auto ja = std::optional<unmasked_fom_jac_action_result_type*>(&unMaskedFomJacAction_);
      fomSystem_.get().residualAndJacobianAction(fomState_, unMaskedFomResidual_, phi, ja);
    }else{
      fomSystem_.get().residualAndJacobianAction(fomState_, unMaskedFomResidual_, phi, {});
    }

    // do masking
    masker_(unMaskedFomResidual_, lspgResidual);
    if (lspgJacobian){
      masker_(unMaskedFomJacAction_, *lspgJacobian.value());
    }

    scaler_(fomState_, lspgResidual, lspgJacobian);
  }

protected:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const MaskerType> masker_;
  mutable unmasked_fom_residual_type unMaskedFomResidual_;
  mutable unmasked_fom_jac_action_result_type unMaskedFomJacAction_;
  PossiblyRefWrapperOperatorScalerType scaler_;
};

}}} // end pressio::rom::impl
#endif  // PRESSIO_ROM_IMPL_LSPG_STEADY_SYSTEM_MASKED_HPP_
