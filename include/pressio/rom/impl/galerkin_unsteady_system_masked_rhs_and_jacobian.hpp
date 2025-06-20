
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_AND_JACOBIAN_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_AND_JACOBIAN_HPP_

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
#endif  // PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_AND_JACOBIAN_HPP_
