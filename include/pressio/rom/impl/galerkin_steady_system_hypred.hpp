
#ifndef PRESSIO_ROM_IMPL_GALERKIN_STEADY_SYSTEM_HYPRED_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_STEADY_SYSTEM_HYPRED_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
hypred Galerkin problem represents:

   hypredOp fom_r(phi x) = 0

- fom_r is the fom "residual"
- phi is the basis

From this we get a "reduced" residual/jacobian:
R = phi^T hypredOp fom_r(phi x)
J = phi^T hypredOp dfom_r/dx(phi x) phi
*/
template <
  class ReducedStateType,
  class ReducedResidualType,
  class ReducedJacobianType,
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReductionOperator
  >
class GalerkinSteadyHypRedSystem
{

  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;

  // deduce from the fom object the type of result of
  // applying the Jacobian to the basis
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createResultOfJacobianActionOn
	     (std::declval<basis_matrix_type const &>()) );

public:
  // required aliases
  using state_type    = ReducedStateType;
  using residual_type = ReducedResidualType;
  using jacobian_type = ReducedJacobianType;

  GalerkinSteadyHypRedSystem() = delete;

  GalerkinSteadyHypRedSystem(const TrialSubspaceType & trialSubspace,
			     const FomSystemType & fomSystem,
			     const HyperReductionOperator & hyperReducer)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      hyperReducer_(hyperReducer),
      fomResidual_(fomSystem.createResidual()),
      fomJacAction_(fomSystem.createResultOfJacobianActionOn(trialSubspace_.get().basisOfTranslatedSpace()))
  {}

  GalerkinSteadyHypRedSystem(GalerkinSteadyHypRedSystem const &) = delete;
  GalerkinSteadyHypRedSystem& operator=(GalerkinSteadyHypRedSystem const&) = delete;
  GalerkinSteadyHypRedSystem(GalerkinSteadyHypRedSystem &&) = default;
  GalerkinSteadyHypRedSystem& operator=(GalerkinSteadyHypRedSystem &&) = default;

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  residual_type createResidual() const{
    // here we assume that the action of hyOp does not
    // produce a reduced residual different than the number of basis.
    // to be precise, we should compute: R = MJOP f(phi x)
    return impl::CreateGalerkinRhs<residual_type>()(trialSubspace_.get().dimension());
  }

  jacobian_type createJacobian() const{
    // here we assume that the reduced jacobian is square matrix
    // defined by num of modes in basis.
    // to be precise, we should compute: J = MJOP df/dx(phi x) phi
    return impl::CreateGalerkinJacobian<jacobian_type>()(trialSubspace_.get().dimension());
  }

  void residualAndJacobian(const state_type & reducedState,
			   residual_type & reducedResidual,
			   std::optional<jacobian_type*> reducedJacobian) const
  {
    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    std::optional<fom_jac_action_result_type *> fomJacActionOpt;
    if (reducedJacobian) {
      fomJacActionOpt = &fomJacAction_;
    }
    fomSystem_.get().residualAndJacobianAction(fomState_, fomResidual_, phi, fomJacActionOpt);

    hyperReducer_(fomResidual_, reducedResidual);
    if (reducedJacobian){
      hyperReducer_(fomJacAction_, *reducedJacobian.value());
    }
  }

  void residual(const state_type & reducedState,
		residual_type & reducedResidual) const
  {
    throw std::runtime_error("missing impl");
  }

  template<class OperandT, class ResultT>
  void applyJacobian(const state_type & reducedState,
		     OperandT const & reducedOperand,
		     ResultT & out) const
  {
    throw std::runtime_error("missing impl");
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReductionOperator> hyperReducer_;
  mutable typename FomSystemType::residual_type fomResidual_;
  mutable fom_jac_action_result_type fomJacAction_;
};

}}}
#endif  // PRESSIO_ROM_IMPL_GALERKIN_STEADY_SYSTEM_HYPRED_HPP_
