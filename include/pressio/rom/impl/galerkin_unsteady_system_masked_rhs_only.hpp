
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_ONLY_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_ONLY_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  masked explicit galerkin system represents:

     d hat{y}/dt = HrOp mask( fom_rhs(phi*hat{y}, ...) )

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- HrOp is the hyper-red operator
*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType
  >
class GalerkinMaskedOdeSystemOnlyRhs
{
  // deduce types
  using unmasked_fom_rhs_type = typename FomSystemType::rhs_type;
  using masked_fom_rhs_type =
    decltype(std::declval<MaskerType const>().createResultOfMaskActionOn
	     (std::declval<unmasked_fom_rhs_type const &>()));

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using rhs_type      = ReducedRhsType;

  GalerkinMaskedOdeSystemOnlyRhs(const TrialSubspaceType & trialSubspace,
				 const FomSystemType & fomSystem,
				 const MaskerType & masker,
				 const HyperReducerType & hyperReducer)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      hyperReducer_(hyperReducer),
      masker_(masker),
      unMaskedFomRhs_(fomSystem.createRhs()),
      maskedFomRhs_(masker.createResultOfMaskActionOn(unMaskedFomRhs_))
  {}


  GalerkinMaskedOdeSystemOnlyRhs(GalerkinMaskedOdeSystemOnlyRhs const &) = delete;
  GalerkinMaskedOdeSystemOnlyRhs& operator=(GalerkinMaskedOdeSystemOnlyRhs const&) = delete;
  GalerkinMaskedOdeSystemOnlyRhs(GalerkinMaskedOdeSystemOnlyRhs &&) = default;
  GalerkinMaskedOdeSystemOnlyRhs& operator=(GalerkinMaskedOdeSystemOnlyRhs &&) = default;

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  rhs_type createRhs() const{
    return impl::CreateGalerkinRhs<rhs_type>()(trialSubspace_.get().dimension());
  }

  void rhs(const state_type & reducedState,
	   const IndVarType & rhsEvaluationTime,
	   rhs_type & reducedRhs) const
  {

    // reconstruct fom state fomState = phi*reducedState
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);
    // evaluate fomRhs
    fomSystem_.get().rhs(fomState_, rhsEvaluationTime, unMaskedFomRhs_);
    // mask it
    masker_(unMaskedFomRhs_, maskedFomRhs_);
    // evaluate reduced rhs
    hyperReducer_(maskedFomRhs_, rhsEvaluationTime, reducedRhs);
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReducerType> hyperReducer_;
  std::reference_wrapper<const MaskerType> masker_;
  mutable unmasked_fom_rhs_type unMaskedFomRhs_;
  mutable masked_fom_rhs_type maskedFomRhs_;
};

}}} // end pressio::rom::impl
#endif  // PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_ONLY_HPP_
