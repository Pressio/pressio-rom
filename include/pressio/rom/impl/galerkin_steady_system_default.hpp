
#ifndef PRESSIO_ROM_IMPL_GALERKIN_STEADY_SYSTEM_DEFAULT_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_STEADY_SYSTEM_DEFAULT_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
default Galerkin problem represents:

  phi^T fom_r(phi x) = 0

- fom_r is the fom "residual"
- phi is the basis

From this we get a "reduced" residual/jacobian:
  R = phi^T fom_r(phi x)
  J = phi^T dfom_r/dx(phi x) phi

*/
template <
  class ReducedStateType,
  class ReducedResidualType,
  class ReducedJacobianType,
  class TrialSubspaceType,
  class FomSystemType
  >
class GalerkinSteadyDefaultSystem
{

  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;

  // Type resulting from applying the FOM Jacobian to the basis
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createResultOfJacobianActionOn
	     (std::declval<basis_matrix_type const &>()) );

public:
  // aliases required by the pressio solvers
  using state_type    = ReducedStateType;
  using residual_type = ReducedResidualType;
  using jacobian_type = ReducedJacobianType;

  GalerkinSteadyDefaultSystem() = delete;

  GalerkinSteadyDefaultSystem(const TrialSubspaceType & trialSubspace,
			      const FomSystemType & fomSystem)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      fomResidual_(fomSystem.createResidual()),
      fomJacAction_(fomSystem.createResultOfJacobianActionOn(trialSubspace_.get().basisOfTranslatedSpace())),
      auxFomVec_(pressio::ops::clone(fomState_)),
      auxFomVec2_(pressio::ops::clone(fomState_))
  {}

  GalerkinSteadyDefaultSystem(GalerkinSteadyDefaultSystem const &) = delete;
  GalerkinSteadyDefaultSystem& operator=(GalerkinSteadyDefaultSystem const&) = delete;
  GalerkinSteadyDefaultSystem(GalerkinSteadyDefaultSystem &&) = default;
  GalerkinSteadyDefaultSystem& operator=(GalerkinSteadyDefaultSystem &&) = default;

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

  void residual(const state_type & reducedState,
		residual_type & reducedResidual) const
  {
    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    fomSystem_.get().residual(fomState_, fomResidual_);
    ::pressio::ops::product(::pressio::transpose(),
			    1, phi, fomResidual_, 0, reducedResidual);
  }

  /*
   * Computes the reduced-order residual and, optionally, the reduced Jacobian
   *
   * Parameters:
   *   - reducedState:       current ROM state vector
   *   - reducedResidual:    output reduced residual vector (to be computed)
   *   - reducedJacobian:    optional output reduced Jacobian matrix
   *
   * This function:
   *   1. Maps the reduced state to the full-order space.
   *   2. Computes the FOM residual and (optionally) the action of the FOM Jacobian on the basis.
   *   3. Projects the FOM residual onto the trial subspace to compute the ROM residual.
   *   4. If requested, projects the FOM Jacobian action result to compute the ROM Jacobian.
   */
  void residualAndJacobian(const state_type & reducedState,
			   residual_type & reducedResidual,
			   std::optional<jacobian_type*> reducedJacobian) const
  {

    // Get a reference to the basis matrix from the trial subspace
    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    // Map the reduced state to the full-order state:
    // fomState = phi * reducedState + reference state
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    // Define scalar types and scaling constants for projections
    using phi_scalar_t = typename ::pressio::Traits<basis_matrix_type>::scalar_type;
    constexpr auto alpha = static_cast<phi_scalar_t>(1);
    using R_scalar_t = typename ::pressio::Traits<residual_type>::scalar_type;
    constexpr auto beta = static_cast<R_scalar_t>(0);

    // Prepare an optional pointer to store the result of Jacobian action, if needed
    std::optional<fom_jac_action_result_type *> fomJacActionOpt;
    if (reducedJacobian) {
      fomJacActionOpt = &fomJacAction_;
    }
    // Compute the FOM residual and optionally the Jacobian applied to phi
    fomSystem_.get().residualAndJacobianAction(fomState_, fomResidual_, phi, fomJacActionOpt);

    // Project the FOM residual to the reduced space: reducedResidual = phi^T * fomResidual
    ::pressio::ops::product(::pressio::transpose(),
			    alpha, phi, fomResidual_, beta, reducedResidual);

    // If the reduced Jacobian is requested
    // compute reducedJacobian = phi^T * (J_fom * phi)
    if (reducedJacobian){
      using J_scalar_t = typename ::pressio::Traits<jacobian_type>::scalar_type;
      constexpr auto beta = static_cast<J_scalar_t>(0);
      ::pressio::ops::product(::pressio::transpose(),
			      ::pressio::nontranspose(),
			      alpha, phi, fomJacAction_, beta,
			      *reducedJacobian.value());
    }
  }

  // these is here for matrix-free methods
  template<class OperandT, class ResultT>
  void applyJacobian(const state_type & reducedState,
		     OperandT const & reducedOperand,
		     ResultT & out) const
  {
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    pressio::ops::set_zero(auxFomVec_);
    trialSubspace_.get().mapFromReducedStateWithoutTranslation(reducedOperand, auxFomVec_);
    pressio::ops::set_zero(auxFomVec2_);
    fomSystem_.get().jacobianAction(fomState_, auxFomVec_, auxFomVec2_);

    using scalar_t = typename ::pressio::Traits<basis_matrix_type>::scalar_type;
    scalar_t alpha{1};
    scalar_t beta{0};
    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    ::pressio::ops::product(::pressio::transpose(),
          alpha, phi, auxFomVec2_, beta, out);
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;

  // FOM state, residual and jacobian action
  mutable typename FomSystemType::state_type fomState_;
  mutable typename FomSystemType::residual_type fomResidual_;
  mutable fom_jac_action_result_type fomJacAction_;

  // Extra auxiliary vectors
  mutable typename FomSystemType::state_type auxFomVec_;
  mutable typename FomSystemType::state_type auxFomVec2_;
};

}}} // end pressio::rom::impl
#endif  // PRESSIO_ROM_IMPL_GALERKIN_STEADY_SYSTEM_DEFAULT_HPP_
