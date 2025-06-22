
#ifndef PRESSIO_ROM_GALERKIN_UNSTEADY_EXPLICIT_HPP_
#define PRESSIO_ROM_GALERKIN_UNSTEADY_EXPLICIT_HPP_

#include "./reduced_operators_traits.hpp"
#include "impl/galerkin_helpers.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_only.hpp"
#include "impl/galerkin_unsteady_system_hypred_rhs_only.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_with_mass_matrix.hpp"
#include "impl/galerkin_unsteady_system_masked_rhs_only.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

/*
 These function overloads construct Galerkin ROM representations for
 unsteady problems.

 Overloads:
  (1) Default Galerkin with and without mass matrix
  (2) Hyper-reduced Galerkin (no mass matrix):
  (3) Masked Galerkin with hyper-reduction:

 Requirements:
   - The trial subspace must satisfy PossiblyAffineTrialColumnSubspace.
   - The FOM system must satisfy different concepts depending on the overload

 Return Type:
  Each function returns an object modeling the Pressio explicit stepper concept.
  The concrete type of the returned object is not part of the public API, but it
  satisfies the requirements of an explicit stepper as defined in pressio/ode.

  The returned stepper can then be used to advance the reduced-order model in time,
  using, for example, the `advance_n_steps` or `advance_to_target_time` methods
  provided by the Pressio ODE integration module.
*/

// -------------------------------------------------------------
// default
// -------------------------------------------------------------

template<
 class TrialSubspaceType, class FomSystemType
 >
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,  /*(1)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem)
{

  static_assert(PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value);
  static_assert(RealValuedSemiDiscreteFom<FomSystemType>::value);
  static_assert(std::is_same<typename TrialSubspaceType::
		full_state_type, typename FomSystemType::state_type>::value);

  impl::valid_scheme_for_explicit_galerkin_else_throw(schemeName, "galerkin_default_explicit");
  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using reduced_rhs_type = impl::explicit_galerkin_default_reduced_rhs_t<TrialSubspaceType>;

  constexpr bool with_mass_matrix = RealValuedSemiDiscreteFomWithMassMatrixAction<
    FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value;
  if constexpr (with_mass_matrix){
    using reduced_mm_type = impl::explicit_galerkin_default_reduced_mass_matrix_t<TrialSubspaceType>;
    using galerkin_system = impl::GalerkinDefaultOdeSystemOnlyRhsAndMassMatrix<
      ind_var_type, reduced_state_type, reduced_rhs_type,
      reduced_mm_type, TrialSubspaceType, FomSystemType>;
    auto gs = std::make_unique<galerkin_system>(trialSpace, fomSystem);
    return ::pressio::ode::create_explicit_stepper<galerkin_system>(schemeName, std::move(gs));
  }
  else{
    using galerkin_system = impl::GalerkinDefaultOdeSystemOnlyRhs<
      ind_var_type, reduced_state_type, reduced_rhs_type, TrialSubspaceType, FomSystemType>;
    auto gs = std::make_unique<galerkin_system>(trialSpace, fomSystem);
    return ::pressio::ode::create_explicit_stepper<galerkin_system>(schemeName, std::move(gs));
  }

}


// -------------------------------------------------------------
// hyper-reduced
// -------------------------------------------------------------

template<
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType>
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,  /*(2)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const HyperReducerType & hyperReducer)
{
  static_assert(PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value);
  static_assert(RealValuedSemiDiscreteFom<FomSystemType>::value);

  impl::valid_scheme_for_explicit_galerkin_else_throw(schemeName, "galerkin_hypred_explicit");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using reduced_rhs_type = impl::explicit_galerkin_default_reduced_rhs_t<TrialSubspaceType>;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinHyperReducedOdeSystemOnlyRhs<
    ind_var_type, reduced_state_type, reduced_rhs_type,
    TrialSubspaceType, FomSystemType, HyperReducerType>;

  auto gs = std::make_unique<galerkin_system>(trialSpace, fomSystem, hyperReducer);
  return ::pressio::ode::create_explicit_stepper<galerkin_system>(schemeName, std::move(gs));
}

// -------------------------------------------------------------
// masked
// -------------------------------------------------------------

template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType>
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,  /*(3)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const MaskerType & masker,
				      const HyperReducerType & hyperReducer)
{
  static_assert(PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value);
  static_assert(RealValuedSemiDiscreteFom<FomSystemType>::value);

  impl::valid_scheme_for_explicit_galerkin_else_throw(schemeName, "galerkin_masked_explicit");
  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using reduced_rhs_type = reduced_state_type;
  // the "system" implements the math
  using galerkin_system = impl::GalerkinMaskedOdeSystemOnlyRhs<
    ind_var_type, reduced_state_type, reduced_rhs_type,
    TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>;

  auto gs = std::make_unique<galerkin_system>(trialSpace, fomSystem, masker, hyperReducer);
  return ::pressio::ode::create_explicit_stepper<galerkin_system>(schemeName, std::move(gs));
}

}}} // end pressio::rom::galerkin
#endif  // PRESSIO_ROM_GALERKIN_UNSTEADY_EXPLICIT_HPP_
