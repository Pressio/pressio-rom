
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_HYPRED_RHS_AND_JACOBIAN_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_HYPRED_RHS_AND_JACOBIAN_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  hypred implicit galerkin system represents:

     d hat{y}/dt = hyperReducer fom_rhs(phi*hat{y}, ...)

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- hyperReducer is the hypred operator
so that it boils down to:

rhs = hyperReducer fom_rhs(phi*hat{y}, ...)
rhs_jacobian = hyperReducer d(fom_rhs(phi*hat{y}, ...))/dy phi

*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class ReducedJacobianType,
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType
  >
class GalerkinHypRedOdeSystemRhsAndJacobian
{

  // deduce from the fom object the type of result of
  // applying the Jacobian to the basis
  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createResultOfJacobianActionOn
	     (std::declval<basis_matrix_type const &>()) );

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using rhs_type      = ReducedRhsType;
  using jacobian_type = ReducedJacobianType;

  GalerkinHypRedOdeSystemRhsAndJacobian(const TrialSubspaceType & trialSubspace,
					const FomSystemType & fomSystem,
					const HyperReducerType & hyperReducer)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      hyperReducer_(hyperReducer),
      fomRhs_(fomSystem.createRhs()),
      fomJacAction_(fomSystem.createResultOfJacobianActionOn(trialSubspace_.get().basisOfTranslatedSpace()))
  {}

  GalerkinHypRedOdeSystemRhsAndJacobian(GalerkinHypRedOdeSystemRhsAndJacobian const &) = delete;
  GalerkinHypRedOdeSystemRhsAndJacobian& operator=(GalerkinHypRedOdeSystemRhsAndJacobian const&) = delete;
  GalerkinHypRedOdeSystemRhsAndJacobian(GalerkinHypRedOdeSystemRhsAndJacobian &&) = default;
  GalerkinHypRedOdeSystemRhsAndJacobian& operator=(GalerkinHypRedOdeSystemRhsAndJacobian &&) = default;

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
    fomSystem_.get().rhs(fomState_, rhsEvaluationTime, fomRhs_);
    hyperReducer_(fomRhs_, rhsEvaluationTime, reducedRhs);
    if (reducedJacobian){
      const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
      fomSystem_.get().applyJacobian(fomState_, phi, rhsEvaluationTime, fomJacAction_);
      hyperReducer_(fomJacAction_, rhsEvaluationTime, *reducedJacobian.value());
    }
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReducerType> hyperReducer_;
  mutable typename FomSystemType::rhs_type fomRhs_;
  mutable fom_jac_action_result_type fomJacAction_;
};

}}} // end pressio::rom::impl
#endif  // PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_HYPRED_RHS_AND_JACOBIAN_HPP_
