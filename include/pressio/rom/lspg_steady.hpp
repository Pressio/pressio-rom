/*
//@HEADER
// ************************************************************************
//
// lspg_steady.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef PRESSIOROM_ROM_LSPG_STEADY_HPP_
#define PRESSIOROM_ROM_LSPG_STEADY_HPP_

#include "./impl/lspg_helpers.hpp"
#include "./impl/lspg_steady_system_default.hpp"
#include "./impl/lspg_steady_system_masked.hpp"
#include "./impl/noop.hpp"

namespace pressio{ namespace rom{ namespace lspg{

/*
 * These overloads construct LSPG ROM representations for steady problems.
 *
 * Each overload returns an object modeling a nonlinear steady system in
 * reduced space that exposes a residual and applyJacobian interface,
 * compatible with Pressio's nonlinear solvers.
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

  // check return type meets nonlinearsolvers::NonlinearSystemFusingResidualAndJacobian
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
  // check return type meets nonlinearsolvers::NonlinearSystemFusingResidualAndJacobian
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
  // check return type meets nonlinearsolvers::NonlinearSystemFusingResidualAndJacobian
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
  // check return type meets nonlinearsolvers::NonlinearSystemFusingResidualAndJacobian
  impl::steady_lspg_static_check_api_return_type<return_type>();

  return return_type(trialSpace, fomSystem, masker, scaler);
}

} // end experimental

}}} // end pressio::rom::lspg
#endif  // PRESSIOROM_ROM_LSPG_STEADY_HPP_
