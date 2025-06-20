
#ifndef PRESSIO_ROM_LSPG_STEADY_HPP_
#define PRESSIO_ROM_LSPG_STEADY_HPP_

#include "./impl/lspg_helpers.hpp"
#include "./impl/lspg_steady_system_default.hpp"
#include "./impl/lspg_steady_system_masked.hpp"
#include "./impl/noop.hpp"

namespace pressio{ namespace rom{ namespace lspg{

/*
 * These overloads construct steady LSPG ROM representations
 * for steady problems.
 *
 * The returned object models a nonlinear steady system in reduced space
 * that exposes a residual and applyJacobian interface, compatible with
 * Pressio's nonlinear solvers (e.g., Newton solvers).
 *
 * Overloads:
 *
 * (1) Default or hyper-reduced LSPG:
 *     - No masking or scaling applied.
 *     - Returns a LspgSteadyDefaultSystem.
 *
 * (2) Masked LSPG:
 *     - Applies a residual masker (e.g., for gappy POD or GNAT).
 *     - Returns a LspgSteadyMaskedSystem.
 *
 * (3) Experimental: Default/hyper-reduced LSPG with scaling:
 *     - Applies a user-defined scaling operator (e.g., weighting).
 *     - Returns a LspgSteadyDefaultSystem with scaling.
 *
 * (4) Experimental: Masked LSPG with scaling:
 *     - Combines masking and residual scaling.
 *     - Returns a LspgSteadyMaskedSystem with scaling.
 *
 * Notes:
 *   - The distinction between "default" and "hyper-reduced" LSPG is internal;
 *     both share the same implementation as long as the projection operator
 *     includes hyper-reduction.
 *   - The returned system types are internal implementation details and should
 *     not be relied on directly. Only the nonlinear system interface they expose
 *     (residual, applyJacobian) is guaranteed.
 *     Indeed, the returned object's type satisifies the
 *     NonlinearSystemFusingResidualAndJacobian concept from the pressio nonlinear solvers.
 *
 * Requirements:
 *   - Trial subspace must satisfy PossiblyAffineTrialColumnSubspace.
 *   - FOM system must satisfy the SteadyFomWithJacobianAction concept.
 */

// -------------------------------------------------------------
// default or hyp-red
// impl-wise, default and hyp-red LSPG are the same
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType>
auto create_steady_problem(const TrialSubspaceType & trialSpace,   /*(1)*/
			   const FomSystemType & fomSystem)

{
  impl::lspg_static_check_trial_and_system(trialSpace, fomSystem);

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = impl::NoOperation<void>;
  using return_type = impl::LspgSteadyDefaultSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, scaler_type>;
  impl::steady_lspg_static_check_api_return_type<return_type>();

  return return_type(trialSpace, fomSystem, scaler_type{});
}

// -------------------------------------------------------------
// masked
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType>
auto create_steady_problem(const TrialSubspaceType & trialSpace,  /*(2)*/
			   const FomSystemType & fomSystem,
			   const MaskerType & masker)
{
  impl::lspg_static_check_trial_and_system(trialSpace, fomSystem);

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = impl::NoOperation<void>;
  using return_type = impl::LspgSteadyMaskedSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, MaskerType, scaler_type>;
  impl::steady_lspg_static_check_api_return_type<return_type>();

  return return_type(trialSpace, fomSystem, masker, scaler_type{});
}


namespace experimental{
// -------------------------------------------------------------
// default or hyp-red with scaling
// impl-wise, default and hyp-red LSPG are the same
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class ScalingOperatorType>
auto create_steady_problem(const TrialSubspaceType & trialSpace,  /*(3)*/
			   const FomSystemType & fomSystem,
			   const ScalingOperatorType & scaler)
{
  impl::lspg_static_check_trial_and_system(trialSpace, fomSystem);

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = std::reference_wrapper<const ScalingOperatorType>;
  using return_type = impl::LspgSteadyDefaultSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, scaler_type>;
  impl::steady_lspg_static_check_api_return_type<return_type>();

  return return_type(trialSpace, fomSystem, scaler);
}

// -------------------------------------------------------------
// masked with scaling
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class ScalingOperatorType>
auto create_steady_problem(const TrialSubspaceType & trialSpace,  /*(4)*/
			   const FomSystemType & fomSystem,
			   const MaskerType & masker,
			   const ScalingOperatorType & scaler)
{
  impl::lspg_static_check_trial_and_system(trialSpace, fomSystem);

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = std::reference_wrapper<const ScalingOperatorType>;
  using return_type = impl::LspgSteadyMaskedSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, MaskerType, scaler_type>;
  impl::steady_lspg_static_check_api_return_type<return_type>();

  return return_type(trialSpace, fomSystem, masker, scaler);
}

} // end experimental

}}} // end pressio::rom::lspg
#endif  // PRESSIO_ROM_LSPG_STEADY_HPP_
