
#ifndef PRESSIO_ROM_GALERKIN_STEADY_HPP_
#define PRESSIO_ROM_GALERKIN_STEADY_HPP_

#include "./impl/galerkin_helpers.hpp"
#include "./impl/galerkin_steady_system_default.hpp"
#include "./impl/galerkin_steady_system_hypred.hpp"
#include "./impl/galerkin_steady_system_masked.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

/*
 * These function overloads construct Galerkin ROM representations for
 * steady problems.
 * Each overload returns an object modeling a nonlinear steady system
 * in the reduced space, which exposes a residual and applyJacobian interface
 * compatible with Pressio's nonlinear solvers.
 *
 * Overloads:
 *
 * (1) Default Galerkin problem
 *
 * (2) Hyper-reduced Galerkin problem
 *
 * (3) Masked-based hyper-reduced Galerkin problem
 *
 * Requirements:
 *   - The trial subspace must satisfy PossiblyAffineTrialColumnSubspace.
 *   - The FOM system must satisfy the SteadyFomWithJacobianAction concept
 *
 * Return Type:
 *   The return type of each overload is an internal implementation detail
 *   and is not part of the public API. It is intentionally opaque and subject to change.
 *   Users should rely only on the fact that the returned object satisfies the interface
 *   of a nonlinear system (i.e., exposes residual and applyJacobian methods, etc).
 *
 * These returned objects can be used directly to compute the solution of
 * the steady ROM, e.g., using the Pressio solvers such as create_newton_solver.
 */


// ------------------------------------------------------------------------
// default
// ------------------------------------------------------------------------

template<
  class TrialSubspaceType,
  class FomSystemType>
auto create_steady_problem(const TrialSubspaceType & trialSpace,   /*(1)*/
			   const FomSystemType & fomSystem)
{

  static_assert(PossiblyAffineTrialColumnSubspace<TrialSubspaceType>::value,
     "You are trying to create a steady galerkin problem but the \
trialSpace does not meet the required PossiblyAffineTrialColumnSubspace concept.");

  static_assert
    (SteadyFomWithJacobianAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
     "You are trying to create a steady galerkin problem but the \
FOM system does not meet the required SteadyFomWithJacobianAction concept.");

  // figure out what types we need to use for the "reduced" system.
  // deduce this from the reduced state the user set for the trial subspace.
  using reduced_state_t = typename TrialSubspaceType::reduced_state_type;
  using reduced_residual_t = impl::steady_galerkin_default_reduced_residual_t<reduced_state_t>;
  using reduced_jacobian_t = impl::steady_galerkin_default_reduced_jacobian_t<reduced_state_t>;
  using return_type = impl::GalerkinSteadyDefaultSystem<
    reduced_state_t, reduced_residual_t, reduced_jacobian_t,
    TrialSubspaceType, FomSystemType>;
  return return_type(trialSpace, fomSystem);
}

// ------------------------------------------------------------------------
// hyper-reduced
// ------------------------------------------------------------------------

template<
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType>
auto create_steady_problem(const TrialSubspaceType & trialSpace,   /*(2)*/
			   const FomSystemType & fomSystem,
			   const HyperReducerType & hyperReducer)
{

  static_assert(PossiblyAffineTrialColumnSubspace<TrialSubspaceType>::value,
     "You are trying to create a steady galerkin problem but the \
trialSpace does not meet the required PossiblyAffineTrialColumnSubspace concept.");

  static_assert
    (SteadyFomWithJacobianAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
     "You are trying to create a steady galerkin problem but the \
FOM system does not meet the required SteadyFomWithJacobianAction concept.");

  using reduced_state_t = typename TrialSubspaceType::reduced_state_type;
  using reduced_residual_t = impl::steady_galerkin_default_reduced_residual_t<reduced_state_t>;
  using reduced_jacobian_t = impl::steady_galerkin_default_reduced_jacobian_t<reduced_state_t>;
  using return_type = impl::GalerkinSteadyHypRedSystem<
    reduced_state_t, reduced_residual_t, reduced_jacobian_t,
    TrialSubspaceType, FomSystemType, HyperReducerType>;
  return return_type(trialSpace, fomSystem, hyperReducer);
}

// ------------------------------------------------------------------------
// masked-based hyper-reduced
// ------------------------------------------------------------------------

template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType>
auto create_steady_problem(const TrialSubspaceType & trialSpace,   /*(3)*/
			   const FomSystemType & fomSystem,
			   const MaskerType & masker,
			   const HyperReducerType & hyperReducer)
{
  static_assert(PossiblyAffineTrialColumnSubspace<TrialSubspaceType>::value,
     "You are trying to create a steady galerkin problem but the \
trialSpace does not meet the required PossiblyAffineTrialColumnSubspace concept.");

  static_assert
    (SteadyFomWithJacobianAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
     "You are trying to create a steady galerkin problem but the \
FOM system does not meet the required SteadyFomWithJacobianAction concept.");

  using reduced_state_t = typename TrialSubspaceType::reduced_state_type;
  using reduced_residual_t = impl::steady_galerkin_default_reduced_residual_t<reduced_state_t>;
  using reduced_jacobian_t = impl::steady_galerkin_default_reduced_jacobian_t<reduced_state_t>;
  using return_type = impl::GalerkinSteadyMaskedSystem<
    reduced_state_t, reduced_residual_t, reduced_jacobian_t,
    TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>;
  return return_type(trialSpace, fomSystem, masker, hyperReducer);
}

}}} // end pressio::rom::galerkin
#endif  // PRESSIO_ROM_GALERKIN_STEADY_HPP_
