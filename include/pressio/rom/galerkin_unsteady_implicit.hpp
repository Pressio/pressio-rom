
#ifndef PRESSIO_ROM_GALERKIN_UNSTEADY_IMPLICIT_HPP_
#define PRESSIO_ROM_GALERKIN_UNSTEADY_IMPLICIT_HPP_

#include "impl/galerkin_helpers.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_and_jacobian_and_mm.hpp"
#include "impl/galerkin_unsteady_system_hypred_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_masked_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_fully_discrete_fom.hpp"
#include "impl/galerkin_unsteady_system_hypred_fully_discrete_fom.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

// -------------------------------------------------------------
// default
// -------------------------------------------------------------

template<
 class TrialSubspaceType, class FomSystemType,
 std::enable_if_t<
   PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value
   && RealValuedSemiDiscreteFomWithJacobianAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value
   && !RealValuedSemiDiscreteFomWithJacobianAndMassMatrixAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value
   && std::is_same<typename TrialSubspaceType::full_state_type,
		   typename FomSystemType::state_type>::value
   , int > = 0
  >
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,   /*(1)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem)
{

  impl::valid_scheme_for_implicit_galerkin_else_throw(schemeName, "galerkin_default_implicit");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types      = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinDefaultOdeSystemRhsAndJacobian<
    ind_var_type, reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType>;

  auto gs = std::make_unique<galerkin_system>(trialSpace, fomSystem);
  return ::pressio::ode::create_implicit_stepper<galerkin_system>(schemeName, std::move(gs));
}

// -------------------------------------------------------------
// default with mass matrix
// -------------------------------------------------------------

template<
  class TrialSubspaceType, class FomSystemType,
  std::enable_if_t<
    PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value
    && RealValuedSemiDiscreteFomWithJacobianAndMassMatrixAction<
         FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value
    && std::is_same<typename TrialSubspaceType::full_state_type,
		    typename FomSystemType::state_type>::value
    , int > = 0
  >
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,   /*(2)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem)
{

  impl::valid_scheme_for_implicit_galerkin_else_throw(schemeName, "galerkin_default_implicit");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types      = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;
  using reduced_mm_type       = typename default_types::reduced_mass_matrix_type;

  using galerkin_system = impl::GalerkinDefaultOdeSystemRhsJacobianMassMatrix<
    ind_var_type, reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, reduced_mm_type, TrialSubspaceType, FomSystemType>;

  auto gs = std::make_unique<galerkin_system>(trialSpace, fomSystem);
  return ::pressio::ode::create_implicit_stepper(schemeName, std::move(gs));
}


// -------------------------------------------------------------
// hyper-reduced
// -------------------------------------------------------------

template<
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType>
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,   /*(3)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const HyperReducerType & hyperReducer)
{

  impl::valid_scheme_for_implicit_galerkin_else_throw(schemeName, "galerkin_hypred_implicit");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types      = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinHypRedOdeSystemRhsAndJacobian<
    ind_var_type, reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType, HyperReducerType>;

  auto gs = std::make_unique<galerkin_system>(trialSpace, fomSystem, hyperReducer);
  return ::pressio::ode::create_implicit_stepper(schemeName, std::move(gs));
}

// -------------------------------------------------------------
// masked
// -------------------------------------------------------------

template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType>
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,   /*(4)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const MaskerType & masker,
				      const HyperReducerType & hyperReducer)
{

  impl::valid_scheme_for_implicit_galerkin_else_throw(schemeName, "galerkin_masked_implicit");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types      = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinMaskedOdeSystemRhsAndJacobian<
    ind_var_type, reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType,
    MaskerType, HyperReducerType>;

  auto gs = std::make_unique<galerkin_system>(trialSpace, fomSystem, masker, hyperReducer);
  return ::pressio::ode::create_implicit_stepper(schemeName, std::move(gs));
}

// -------------------------------------------------------------
// fully-discrete
// -------------------------------------------------------------

template<
  std::size_t TotalNumberOfDesiredStates,
  class TrialSubspaceType,
  class FomSystemType>
auto create_unsteady_implicit_problem(const TrialSubspaceType & trialSpace,    /*(5)*/
				      const FomSystemType & fomSystem)
{

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type    = typename TrialSubspaceType::reduced_state_type;
  using default_types         = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinDefaultFullyDiscreteSystem<
    TotalNumberOfDesiredStates, ind_var_type,
    reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType>;

  auto gs = std::make_unique<galerkin_system>(trialSpace, fomSystem);
  return ::pressio::ode::create_implicit_stepper<
    TotalNumberOfDesiredStates>(std::move(gs));
}

template<
  std::size_t TotalNumberOfDesiredStates,
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType>
auto create_unsteady_implicit_problem(const TrialSubspaceType & trialSpace,
                                      const FomSystemType & fomSystem,
                                      const HyperReducerType & hyperReducer)
{

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type        = typename TrialSubspaceType::reduced_state_type;
  using default_types             = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinHypRedFullyDiscreteSystem<
    TotalNumberOfDesiredStates, independent_variable_type,
    reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType,
    HyperReducerType>;

  auto gs = std::make_unique<galerkin_system>(trialSpace, fomSystem, hyperReducer);
  return ::pressio::ode::create_implicit_stepper<
    TotalNumberOfDesiredStates>(std::move(gs));
}

}}} // end pressio::rom::galerkin
#endif  // PRESSIO_ROM_GALERKIN_UNSTEADY_IMPLICIT_HPP_
