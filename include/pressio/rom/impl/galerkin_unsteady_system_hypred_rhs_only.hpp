
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_HYPRED_RHS_ONLY_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_HYPRED_RHS_ONLY_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  hyp-red explicit galerkin system represents:

     d hat{y}/dt = hyperReducer fom_rhs(phi*hat{y}, ...)

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- hyperReducer is the hyper-red operator
*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType
  >
class GalerkinHyperReducedOdeSystemOnlyRhs
{

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using rhs_type		  = ReducedRhsType;

  GalerkinHyperReducedOdeSystemOnlyRhs(const TrialSubspaceType & trialSubspace,
				       const FomSystemType & fomSystem,
				       const HyperReducerType & hyperReducer)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      hyperReducer_(hyperReducer),
      fomRhs_(fomSystem.createRhs())
  {}

  GalerkinHyperReducedOdeSystemOnlyRhs(GalerkinHyperReducedOdeSystemOnlyRhs const &) = delete;
  GalerkinHyperReducedOdeSystemOnlyRhs& operator=(GalerkinHyperReducedOdeSystemOnlyRhs const&) = delete;
  GalerkinHyperReducedOdeSystemOnlyRhs(GalerkinHyperReducedOdeSystemOnlyRhs &&) = default;
  GalerkinHyperReducedOdeSystemOnlyRhs& operator=(GalerkinHyperReducedOdeSystemOnlyRhs &&) = default;

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
    fomSystem_.get().rhs(fomState_, rhsEvaluationTime, fomRhs_);
    // evaluate reduced rhs
    hyperReducer_(fomRhs_, rhsEvaluationTime, reducedRhs);
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReducerType> hyperReducer_;
  mutable typename FomSystemType::rhs_type fomRhs_;
};

}}} // end pressio::rom::impl
#endif  // PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_HYPRED_RHS_ONLY_HPP_
