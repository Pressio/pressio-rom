
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_FULLY_DISCRETE_FOM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_FULLY_DISCRETE_FOM_HPP_

#include "./galerkin_unsteady_fom_states_manager.hpp"

namespace pressio{ namespace rom{ namespace impl{

/*
  hat{R} := phi^T fom_R(t_n, y_n+1, t_n, ...)
  hat{J} := phi^T fom_J phi

- fom_R is the fom discrete time residual
- fom_J is the fom discrete time jacobian
*/

template <
  std::size_t n,
  class IndVarType,
  class ReducedStateType,
  class ReducedResidualType,
  class ReducedJacobianType,
  class TrialSubspaceType,
  class FomSystemType
  >
class GalerkinDefaultFullyDiscreteSystem
{
  using raw_step_type = typename ::pressio::ode::StepCount::value_type;
  static_assert(std::is_signed<raw_step_type>::value, "");

  // deduce from the fom object the type of result of
  // applying the Jacobian to the basis
  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createResultOfDiscreteTimeJacobianActionOn
	     (std::declval<basis_matrix_type const &>()) );

  using fom_state_type = typename TrialSubspaceType::full_state_type;

public:
  // required
  using independent_variable_type = IndVarType;
  using state_type	          = ReducedStateType;
  using discrete_residual_type    = ReducedResidualType;
  using discrete_jacobian_type    = ReducedJacobianType;

  GalerkinDefaultFullyDiscreteSystem(const TrialSubspaceType & trialSubspace,
				     const FomSystemType & fomSystem)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomStatesManager_(create_galerkin_fom_states_manager<n>(trialSubspace)),
      fomResidual_(fomSystem.createDiscreteTimeResidual()),
      fomJacAction_(fomSystem.createResultOfDiscreteTimeJacobianActionOn(trialSubspace_.get().basisOfTranslatedSpace()))
  {}

  state_type createState() const{
    // this needs to instantiate the reduced state
    return trialSubspace_.get().createReducedState();
  }

  discrete_residual_type createDiscreteResidual() const{
    // this needs to instantiate the reduced residual
    return impl::CreateGalerkinRhs<discrete_residual_type>()(trialSubspace_.get().dimension());
  }

  discrete_jacobian_type createDiscreteJacobian() const{
    // this needs to instantiate the reduced jacobian
    return impl::CreateGalerkinJacobian<discrete_jacobian_type>()(trialSubspace_.get().dimension());
  }

  template<class StepIntType, std::size_t _n = n>
  std::enable_if_t< (_n == 2)>
  preStepHook(StepIntType stepNumber,
	      IndVarType time, IndVarType dt,
	      const state_type & galerkin_state_np1,
	      const state_type & galerkin_state_n) const
  {
    /* this method is called once before starting a step, so we need to update
       the proposed state and the previous state */
    fomStatesManager_.reconstructAtWithoutStencilUpdate(galerkin_state_np1,
							::pressio::ode::nPlusOne());
    fomStatesManager_.reconstructAtWithoutStencilUpdate(galerkin_state_n,
							::pressio::ode::n());

    static constexpr bool hashook = ::pressio::ode::has_const_pre_step_hook_method<
      mpl::remove_cvref_t<FomSystemType>, _n, StepIntType, IndVarType, fom_state_type
      >::value;
    if constexpr(hashook){
      const auto & ynp1 = fomStatesManager_(::pressio::ode::nPlusOne());
      const auto & yn   = fomStatesManager_(::pressio::ode::n());
      fomSystem_.get().preStepHook(stepNumber, time, dt, ynp1, yn);
    }
  }

  template<typename step_t, std::size_t _n = n>
  std::enable_if_t< (_n==2) >
  discreteResidualAndJacobian(const step_t & currentStepNumber,
			      const independent_variable_type & time_np1,
			      const independent_variable_type & dt,
			      discrete_residual_type & galerkinResidual,
			      std::optional<discrete_jacobian_type *> galerkinJacobian,
			      const state_type & galerkin_state_np1,
			      const state_type & galerkin_state_n) const
  {
    /* the FOM state corresponding to the new predicted state has to be
     * recomputed every time since this method is called at every update
     * of the solution inside a non-linear solve loop
     * Note that here we don't reconstruct at "n" because that state was
     * reconstructed inside the preStepHook and has not changed inside the
     * solver loop
     */
    fomStatesManager_.reconstructAtWithoutStencilUpdate(galerkin_state_np1,
							::pressio::ode::nPlusOne());

    const auto & ynp1 = fomStatesManager_(::pressio::ode::nPlusOne());
    const auto & yn   = fomStatesManager_(::pressio::ode::n());
    const bool computeJacobian = (bool) galerkinJacobian;

    try{
      queryFomOperators(currentStepNumber, time_np1, dt, computeJacobian, ynp1, yn);
    }
    catch (::pressio::eh::DiscreteTimeResidualFailureUnrecoverable const & e){
      throw ::pressio::eh::ResidualEvaluationFailureUnrecoverable();
    }

    computeReducedOperators(galerkinResidual, galerkinJacobian);
  }

private:
  template<typename step_t, class ...States>
  void queryFomOperators(const step_t & currentStepNumber,
			 const independent_variable_type & time_np1,
			 const independent_variable_type & dt,
			 bool computeJacobian,
			 States && ... states) const
  {
    const auto phi = trialSubspace_.get().basisOfTranslatedSpace();

    if (computeJacobian){
      using op_ja_t = std::optional<fom_jac_action_result_type*>;
      fomSystem_.get().discreteTimeResidualAndJacobianAction(currentStepNumber, time_np1,
							     dt, fomResidual_, phi,
							     op_ja_t{&fomJacAction_},
							     std::forward<States>(states)...);
    }
    else{
      fomSystem_.get().discreteTimeResidualAndJacobianAction(currentStepNumber, time_np1,
							     dt, fomResidual_, phi, {},
							     std::forward<States>(states)...);
    }
  }

  void computeReducedOperators(discrete_residual_type & galerkinResidual,
			      std::optional<discrete_jacobian_type *> galerkinJacobian) const
  {
    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    using phi_scalar_t = typename ::pressio::Traits<basis_matrix_type>::scalar_type;
    constexpr auto alpha = static_cast<phi_scalar_t>(1);

    using residual_scalar_t = typename ::pressio::Traits<discrete_residual_type>::scalar_type;
    constexpr auto beta = static_cast<residual_scalar_t>(0);;
    ::pressio::ops::product(::pressio::transpose(),
			    alpha, phi, fomResidual_,
			    beta, galerkinResidual);

    if (galerkinJacobian){
      constexpr auto beta = static_cast<residual_scalar_t>(0);;
      ::pressio::ops::product(::pressio::transpose(),
			      ::pressio::nontranspose(),
			      alpha, phi, fomJacAction_,
			      beta,
			      *galerkinJacobian.value());
    }
  }

protected:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable GalerkinFomStatesManager<TrialSubspaceType> fomStatesManager_;
  mutable typename FomSystemType::discrete_residual_type fomResidual_;
  mutable fom_jac_action_result_type fomJacAction_;
};

}}} // end pressio::rom::impl
#endif  // PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_FULLY_DISCRETE_FOM_HPP_
