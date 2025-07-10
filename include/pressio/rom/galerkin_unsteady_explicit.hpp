/*
//@HEADER
// ************************************************************************
//
// galerkin_unsteady_explicit.hpp
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

#ifndef PRESSIOROM_ROM_GALERKIN_UNSTEADY_EXPLICIT_HPP_
#define PRESSIOROM_ROM_GALERKIN_UNSTEADY_EXPLICIT_HPP_

#include "./reduced_operators_traits.hpp"
#include "impl/galerkin_helpers.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_only.hpp"
#include "impl/galerkin_unsteady_system_hypred_rhs_only.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_with_mass_matrix.hpp"
#include "impl/galerkin_unsteady_system_masked_rhs_only.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

/*
 These function overloads construct Galerkin ROM representations for
 unsteady problems to be integrated explicitly.

 Overloads:
  (1) Default Galerkin with and without mass matrix
  (2) Hyper-reduced Galerkin (no mass matrix):
  (3) Masked Galerkin with hyper-reduction:

 Requirements:
   - The trial subspace must satisfy PossiblyAffineTrialColumnSubspace.
   - The FOM system must satisfy different concepts depending on the overload

 Return Type:
  Each function returns a Pressio explicit stepper.
  The returned stepper can then be used to advance the reduced-order model in time,
  using, for example, the `advance_n_steps` or `advance_to_target_time` methods
  provided by the Pressio ODE integration module.
*/

// -------------------------------------------------------------
// default
// -------------------------------------------------------------

template<
  class TrialSubspaceType,
  class FomSystemType
 >
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,  /*(1)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem)
{

  // check that the trial subspace meets the right concept
  static_assert(PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value);
  // check that the fom system meets the right concept
  static_assert(RealValuedSemiDiscreteFom<FomSystemType>::value);
  // check that the trial space reconstructs a fom state that is compatible
  // with the fom state used in the fom system class
  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value);

  // explicit galerkin requires an explicit ode scheme
  impl::valid_scheme_for_explicit_galerkin_else_throw(schemeName, "galerkin_default_explicit");

  // figure out the type of the independent variable
  using ind_var_type = typename FomSystemType::time_type;

  // figure out what types we need to use for the "reduced" system.
  // deduce this from the reduced state the user set for the trial subspace.
  // for example, if its an Eigen vector, then all the reduced "things"
  // will be set accordingly so that the reduced system is fully represented
  // with Eigen data structures.
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types = ExplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_rhs_type = typename default_types::reduced_rhs_type;

  // if the user-provided system exposes a mass matrix, then we need to return
  // the proper class
  constexpr bool with_mass_matrix = RealValuedSemiDiscreteFomWithMassMatrixAction<
    FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value;
  if constexpr (with_mass_matrix){
    // set the type of the reduced mass matrix
    using reduced_mm_type = typename default_types::reduced_mass_matrix_type;

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
  // check that the trial subspace meets the right concept
  static_assert(PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value);
  // check that the fom system meets the right concept
  static_assert(RealValuedSemiDiscreteFom<FomSystemType>::value);
  // check that the trial space reconstructs a fom state that is compatible
  // with the fom state used in the fom system class
  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value);

  // explicit galerkin requires an explicit ode scheme
  impl::valid_scheme_for_explicit_galerkin_else_throw(schemeName, "galerkin_hypred_explicit");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types = ExplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_rhs_type = typename default_types::reduced_rhs_type;

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
  // the "system" implements the math
  static_assert(PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value);
  static_assert(RealValuedSemiDiscreteFom<FomSystemType>::value);
  // check that the trial space reconstructs a fom state that is compatible
  // with the fom state used in the fom system class
  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value);

  // explicit galerkin requires an explicit ode scheme
  impl::valid_scheme_for_explicit_galerkin_else_throw(schemeName, "galerkin_masked_explicit");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types = ExplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_rhs_type = typename default_types::reduced_rhs_type;

  using galerkin_system = impl::GalerkinMaskedOdeSystemOnlyRhs<
    ind_var_type, reduced_state_type, reduced_rhs_type,
    TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>;

  auto gs = std::make_unique<galerkin_system>(trialSpace, fomSystem, masker, hyperReducer);
  return ::pressio::ode::create_explicit_stepper<galerkin_system>(schemeName, std::move(gs));
}

}}} // end pressio::rom::galerkin
#endif  // PRESSIOROM_ROM_GALERKIN_UNSTEADY_EXPLICIT_HPP_
