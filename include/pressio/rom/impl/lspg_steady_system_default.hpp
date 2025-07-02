
#ifndef PRESSIOROM_ROM_IMPL_LSPG_STEADY_SYSTEM_DEFAULT_HPP_
#define PRESSIOROM_ROM_IMPL_LSPG_STEADY_SYSTEM_DEFAULT_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
LSPG steady default represents:

  min_x ||fom_r(phi x)||

- fom_r is the fom "residual"
- phi is the basis
*/
template <
  class ReducedStateType,
  class TrialSubspaceType,
  class FomSystemType,
  class PossiblyRefWrapperOperatorScalerType
  >
class LspgSteadyDefaultSystem
{

  // need to deduce the type of the action of the fom jacobian
  // which becomes the jacobian_type of the problem
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createResultOfJacobianActionOn
	     (std::declval<typename TrialSubspaceType::basis_matrix_type const &>())
	     );

public:
  // aliases required by the pressio solvers
  using state_type    = ReducedStateType;
  using residual_type = typename FomSystemType::residual_type;
  using jacobian_type = fom_jac_action_result_type;

  // here _RawScalerType must be a template because it is the raw scaler type
  // which can be forwarded to the reference wrapper potentially
  template<class _RawScalerType>
  LspgSteadyDefaultSystem(const TrialSubspaceType & trialSubspace,
			  const FomSystemType & fomSystem,
			  _RawScalerType && scaler)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      scaler_(std::forward<_RawScalerType>(scaler))
  {}

  LspgSteadyDefaultSystem(LspgSteadyDefaultSystem const &) = delete;
  LspgSteadyDefaultSystem& operator=(LspgSteadyDefaultSystem const&) = delete;
  LspgSteadyDefaultSystem(LspgSteadyDefaultSystem &&) = default;
  LspgSteadyDefaultSystem& operator=(LspgSteadyDefaultSystem &&) = default;

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  residual_type createResidual() const{
    return fomSystem_.get().createResidual();
  }

  jacobian_type createJacobian() const{
    return fomSystem_.get().createResultOfJacobianActionOn
      (trialSubspace_.get().basisOfTranslatedSpace());
  }

  void residualAndJacobian(const state_type & lspgState,
			   residual_type & lspgResidual,
			   std::optional<jacobian_type *> lspgJacobian) const
  {
    trialSubspace_.get().mapFromReducedState(lspgState, fomState_);

    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    fomSystem_.get().residualAndJacobianAction(fomState_,
					       lspgResidual,
					       phi, lspgJacobian);
    scaler_(fomState_, lspgResidual, lspgJacobian);
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  PossiblyRefWrapperOperatorScalerType scaler_;
};


}}} // end pressio::rom::impl
#endif  // PRESSIOROM_ROM_IMPL_LSPG_STEADY_SYSTEM_DEFAULT_HPP_
