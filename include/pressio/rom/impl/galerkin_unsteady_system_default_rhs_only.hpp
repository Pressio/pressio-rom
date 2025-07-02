
#ifndef PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_DEFAULT_RHS_ONLY_HPP_
#define PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_DEFAULT_RHS_ONLY_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  default explicit galerkin system represents:

     d hat{y}/dt = phi^T fom_rhs(phi*hat{y}, ...)

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis

*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class TrialSubspaceType,
  class FomSystemType
  >
class GalerkinDefaultOdeSystemOnlyRhs
{
  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using rhs_type		  = ReducedRhsType;

  GalerkinDefaultOdeSystemOnlyRhs(const TrialSubspaceType & trialSubspace,
				  const FomSystemType & fomSystem)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      fomRhs_(fomSystem.createRhs())
  {}

  GalerkinDefaultOdeSystemOnlyRhs(GalerkinDefaultOdeSystemOnlyRhs const &) = delete;
  GalerkinDefaultOdeSystemOnlyRhs& operator=(GalerkinDefaultOdeSystemOnlyRhs const&) = delete;
  GalerkinDefaultOdeSystemOnlyRhs(GalerkinDefaultOdeSystemOnlyRhs &&) = default;
  GalerkinDefaultOdeSystemOnlyRhs& operator=(GalerkinDefaultOdeSystemOnlyRhs &&) = default;

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

    // compute the reduced rhs
    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    using phi_scalar_t = typename ::pressio::Traits<basis_matrix_type>::scalar_type;
    constexpr auto alpha = static_cast<phi_scalar_t>(1);
    using rhs_scalar_t = typename ::pressio::Traits<rhs_type>::scalar_type;
    constexpr auto beta = static_cast<rhs_scalar_t>(0);
    ::pressio::ops::product(::pressio::transpose(),
			    alpha, phi, fomRhs_,
			    beta, reducedRhs);
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  mutable typename FomSystemType::rhs_type fomRhs_;
};

}}} // end pressio::rom::impl
#endif  // PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_DEFAULT_RHS_ONLY_HPP_
