
#ifndef PRESSIO_ROM_LSPG_UNSTEADY_HPP_
#define PRESSIO_ROM_LSPG_UNSTEADY_HPP_

#include "./impl/lspg_helpers.hpp"
#include "./impl/lspg_unsteady_fom_states_manager.hpp"
#include "./impl/lspg_unsteady_rj_policy_default.hpp"
#include "./impl/lspg_unsteady_rj_policy_hypred.hpp"
#include "./impl/lspg_unsteady_fully_discrete_system.hpp"
#include "./impl/lspg_unsteady_mask_decorator.hpp"
#include "./impl/lspg_unsteady_scaling_decorator.hpp"
#include "./impl/lspg_unsteady_problem.hpp"
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./impl/lspg_unsteady_reconstructor.hpp"
#endif

namespace pressio{ namespace rom{ namespace lspg{

// -------------------------------------------------------------
// default
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType>
auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,    /*(1)*/
			     const TrialSubspaceType & trialSpace,
			     const FomSystemType & fomSystem)
{

  impl::valid_scheme_for_lspg_else_throw(schemeName);

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using lspg_residual_type = typename FomSystemType::rhs_type;
  using lspg_jacobian_type =
    decltype(
	     fomSystem.createResultOfJacobianActionOn(trialSpace.basisOfTranslatedSpace())
	     );

  // defining an unsteady lspg problem boils down to
  // defining a "custom" residual and jacobian policy
  // since for lspg we need to customize how we do time stepping
  using rj_policy_type = impl::LspgUnsteadyResidualJacobianPolicy<
    ind_var_type, reduced_state_type,
    lspg_residual_type, lspg_jacobian_type,
    TrialSubspaceType, FomSystemType>;

  using return_type = impl::LspgUnsteadyProblemSemiDiscreteAPI<TrialSubspaceType, rj_policy_type>;
  return return_type(schemeName, trialSpace, fomSystem);
}

// -------------------------------------------------------------
// masked
// -------------------------------------------------------------

template<
  class TrialSubspaceType, class FomSystemType, class MaskerType,
  std::enable_if_t<
    PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value
    && RealValuedSemiDiscreteFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value
    && std::is_same<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>::value
    && MaskableWith<typename FomSystemType::rhs_type, MaskerType>::value
    && MaskableWith<impl::fom_jac_action_on_trial_space_t<FomSystemType, TrialSubspaceType>, MaskerType>::value
    , int> = 0
  >
auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,    /*(2)*/
			     const TrialSubspaceType & trialSpace,
			     const FomSystemType & fomSystem,
			     const MaskerType & masker)
{

  impl::valid_scheme_for_lspg_else_throw(schemeName);

  using ind_var_type   = typename FomSystemType::time_type;
  using reduced_state_type          = typename TrialSubspaceType::reduced_state_type;
  using lspg_unmasked_residual_type = typename FomSystemType::rhs_type;
  using lspg_unmasked_jacobian_type = typename TrialSubspaceType::basis_matrix_type;
  using lspg_residual_type =
    decltype(std::declval<MaskerType const>().createResultOfMaskActionOn
	     (std::declval<lspg_unmasked_residual_type const &>()));

  using lspg_jacobian_type =
    decltype(std::declval<MaskerType const>().createResultOfMaskActionOn
	     (std::declval<lspg_unmasked_jacobian_type const &>()));

  using rj_policy_type =
    impl::LspgMaskDecorator<
      MaskerType, lspg_residual_type, lspg_jacobian_type,
      impl::LspgUnsteadyResidualJacobianPolicy<
	ind_var_type, reduced_state_type,
	lspg_unmasked_residual_type, lspg_unmasked_jacobian_type,
	TrialSubspaceType, FomSystemType
	>
    >;

  using return_type = impl::LspgUnsteadyProblemSemiDiscreteAPI<TrialSubspaceType, rj_policy_type>;
  return return_type(schemeName, trialSpace, fomSystem, masker);
}

// -------------------------------------------------------------
// hyp-red
// -------------------------------------------------------------

template<class TrialSubspaceType, class FomSystemType, class HypRedUpdaterType,
  std::enable_if_t<
    PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value
    && RealValuedSemiDiscreteFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value
    && std::is_same<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>::value
    && !MaskableWith<typename FomSystemType::rhs_type, HypRedUpdaterType>::value
    , int> = 0
  >
auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,    /*(3)*/
			     const TrialSubspaceType & trialSpace,
			     const FomSystemType & fomSystem,
			     const HypRedUpdaterType & hypRedUpdater)
{

  impl::valid_scheme_for_lspg_else_throw(schemeName);

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using lspg_residual_type = typename FomSystemType::rhs_type;
  using lspg_jacobian_type =
    decltype(
	     fomSystem.createResultOfJacobianActionOn(trialSpace.basisOfTranslatedSpace())
	     );

  using rj_policy_type = impl::LspgUnsteadyResidualJacobianPolicyHypRed<
    ind_var_type, reduced_state_type,
    lspg_residual_type, lspg_jacobian_type,
    TrialSubspaceType, FomSystemType, HypRedUpdaterType>;

  using return_type = impl::LspgUnsteadyProblemSemiDiscreteAPI<TrialSubspaceType, rj_policy_type>;
  return return_type(schemeName, trialSpace, fomSystem, hypRedUpdater);
}


namespace experimental{

// -------------------------------------------------------------
// default with scaling
// -------------------------------------------------------------

template<
  class TrialSubspaceType,
  class FomSystemType,
  class ScalingOperatorType>
auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,    /*(4)*/
			     const TrialSubspaceType & trialSpace,
			     const FomSystemType & fomSystem,
			     const ScalingOperatorType & scaler)
{

  impl::valid_scheme_for_lspg_else_throw(schemeName);

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using lspg_residual_type = typename FomSystemType::rhs_type;
  using lspg_jacobian_type =
    decltype(
	     fomSystem.createResultOfJacobianActionOn(trialSpace.basisOfTranslatedSpace())
	     );

  using rj_policy_type =
    impl::LspgScalingDecorator<
      ScalingOperatorType, lspg_residual_type, lspg_jacobian_type, TrialSubspaceType,
      impl::LspgUnsteadyResidualJacobianPolicy<
	ind_var_type, reduced_state_type, lspg_residual_type,
	lspg_jacobian_type, TrialSubspaceType, FomSystemType
	>
    >;

  using return_type = impl::LspgUnsteadyProblemSemiDiscreteAPI<TrialSubspaceType, rj_policy_type>;
  return return_type(schemeName, trialSpace, fomSystem, scaler);
}

// -------------------------------------------------------------
// hyper-reduced with scaling
// -------------------------------------------------------------

template<
  class TrialSubspaceType,
  class FomSystemType,
  class HypRedUpdaterType,
  class ScalingOperatorType>
auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,    /*(5)*/
			     const TrialSubspaceType & trialSpace,
			     const FomSystemType & fomSystem,
			     const HypRedUpdaterType & hypRedUpdater,
			     const ScalingOperatorType & scaler)
{

  impl::valid_scheme_for_lspg_else_throw(schemeName);

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using lspg_residual_type = typename FomSystemType::rhs_type;
  using lspg_jacobian_type =
    decltype(
	     fomSystem.createResultOfJacobianActionOn(trialSpace.basisOfTranslatedSpace())
	     );

  using rj_policy_type =
    impl::LspgScalingDecorator<
      ScalingOperatorType, lspg_residual_type, lspg_jacobian_type, TrialSubspaceType,
      impl::LspgUnsteadyResidualJacobianPolicyHypRed<
	ind_var_type, reduced_state_type, lspg_residual_type,
	lspg_jacobian_type, TrialSubspaceType, FomSystemType, HypRedUpdaterType
	>
    >;

  using return_type = impl::LspgUnsteadyProblemSemiDiscreteAPI<TrialSubspaceType, rj_policy_type>;
  return return_type(schemeName, trialSpace, fomSystem, scaler, hypRedUpdater);
}


} //end namespace experimental


// -------------------------------------------------------------
// fully-discrete
// -------------------------------------------------------------

template<
  std::size_t TotalNumberOfDesiredStates,
  class TrialSubspaceType,
  class FomSystemType>
auto create_unsteady_problem(const TrialSubspaceType & trialSpace,     /*(6)*/
			     const FomSystemType & fomSystem)
{
  static_assert(TotalNumberOfDesiredStates == 2,
		"lspg::create_unsteady_problem currently only supports 2 total states");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using lspg_residual_type = typename FomSystemType::discrete_residual_type;
  using lspg_jacobian_type =
    decltype(
	     fomSystem.createResultOfDiscreteTimeJacobianActionOn
	     (trialSpace.basisOfTranslatedSpace())
	     );

  using system_type = impl::LspgFullyDiscreteSystem<
    TotalNumberOfDesiredStates, ind_var_type, reduced_state_type,
    lspg_residual_type, lspg_jacobian_type,
    TrialSubspaceType, FomSystemType>;

  using return_type = impl::LspgUnsteadyProblemFullyDiscreteAPI<
    TotalNumberOfDesiredStates, TrialSubspaceType, system_type>;
  return return_type(trialSpace, fomSystem);
}


#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<
  class TrialSubspaceType,
  std::enable_if_t<
    PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value,
    int > = 0
  >
auto create_reconstructor(const TrialSubspaceType & trialSpace)
{
  using T = impl::LspgReconstructor<TrialSubspaceType>;
  return T(trialSpace);
}
#endif

}}} // end pressio::rom::lspg
#endif  // PRESSIO_ROM_LSPG_UNSTEADY_HPP_
