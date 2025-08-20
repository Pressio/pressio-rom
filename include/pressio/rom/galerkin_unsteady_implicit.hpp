/*
//@HEADER
// ************************************************************************
//
// galerkin_unsteady_implicit.hpp
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

#ifndef PRESSIOROM_ROM_GALERKIN_UNSTEADY_IMPLICIT_HPP_
#define PRESSIOROM_ROM_GALERKIN_UNSTEADY_IMPLICIT_HPP_

#include "impl/galerkin_helpers.hpp"
#include "impl/galerkin_unsteady_fom_states_manager.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_and_jacobian_and_mm.hpp"
#include "impl/galerkin_unsteady_system_hypred_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_masked_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_fully_discrete_fom.hpp"
#include "impl/galerkin_unsteady_system_hypred_fully_discrete_fom.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

/*
 These function overloads construct Galerkin ROM representations for
 unsteady problems to be integrated implicitly.

 Overloads:
  (1) Default Galerkin with and without mass matrix
  (2) Hyper-reduced Galerkin (no mass matrix):
  (3) Masked Galerkin with hyper-reduction:
  (4) fully discrete systems

 Requirements:
   - The trial subspace must satisfy PossiblyAffineTrialColumnSubspace.
   - The FOM system must satisfy different concepts depending on the overload

 Return Type:
  Each function returns a Pressio implicit stepper.
  The returned stepper can then be used in conjunction with a pressio
  Newton solver to advance the reduced-order model in time,
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
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,   /*(1)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem)
{

  // check that the trial subspace meets the right concept
  static_assert(PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value);
  // check that the fom system meets the right concept
  static_assert(RealValuedSemiDiscreteFomWithJacobianAction<
		FomSystemType,
		typename TrialSubspaceType::basis_matrix_type>::value);
  // check that the trial space reconstructs a fom state that is compatible
  // with the fom state used in the fom system class
  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value);

  // implicit galerkin requires an implicit ode scheme
  impl::valid_scheme_for_implicit_galerkin_else_throw(schemeName, "galerkin_default_implicit");

  // figure out what types we need to use for the "reduced" system.
  // deduce this from the reduced state the user set for the trial subspace.
  // for example, if its an Eigen vector, then all the reduced "things"
  // will be set accordingly so that the reduced system is fully represented
  // with Eigen data structures.
  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // if the user-provided system exposes a mass matrix, then we need to return
  // the proper class
  constexpr bool with_mass_matrix = RealValuedSemiDiscreteFomWithJacobianAndMassMatrixAction<
    FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value;

  if constexpr (with_mass_matrix){
    // set the type of the reduced mass matrix
    using reduced_mm_type       = typename default_types::reduced_mass_matrix_type;

    using galerkin_system = impl::GalerkinDefaultOdeSystemRhsJacobianMassMatrix<
      ind_var_type, reduced_state_type, reduced_residual_type,
      reduced_jacobian_type, reduced_mm_type, TrialSubspaceType, FomSystemType>;
    auto gs = std::make_unique<galerkin_system>(trialSpace, fomSystem);
    return ::pressio::ode::create_implicit_stepper(schemeName, std::move(gs));
  }
  else{
    using galerkin_system = impl::GalerkinDefaultOdeSystemRhsAndJacobian<
      ind_var_type, reduced_state_type, reduced_residual_type,
      reduced_jacobian_type, TrialSubspaceType, FomSystemType>;
    auto gs = std::make_unique<galerkin_system>(trialSpace, fomSystem);
    return ::pressio::ode::create_implicit_stepper(schemeName, std::move(gs));
  }
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

  // check that the trial subspace meets the right concept
  static_assert(PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value);
  // check that the fom system meets the right concept
  static_assert(RealValuedSemiDiscreteFomWithJacobianAction<
		FomSystemType,
		typename TrialSubspaceType::basis_matrix_type>::value);
  // check that the trial space reconstructs a fom state that is compatible
  // with the fom state used in the fom system class
  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value);

  // implicit galerkin requires an implicit ode scheme
  impl::valid_scheme_for_implicit_galerkin_else_throw(schemeName, "galerkin_default_implicit");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

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

  // check that the trial subspace meets the right concept
  static_assert(PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value);
  // check that the fom system meets the right concept
  static_assert(RealValuedSemiDiscreteFomWithJacobianAction<
		FomSystemType,
		typename TrialSubspaceType::basis_matrix_type>::value);
  // check that the trial space reconstructs a fom state that is compatible
  // with the fom state used in the fom system class
  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value);

  // implicit galerkin requires an implicit ode scheme
  impl::valid_scheme_for_implicit_galerkin_else_throw(schemeName, "galerkin_default_implicit");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

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

  // check that the trial subspace meets the right concept
  static_assert(PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value);
  // check that the fom system meets the right concept
  static_assert(FullyDiscreteSystemWithJacobianAction<
		FomSystemType,
		TotalNumberOfDesiredStates,
		typename TrialSubspaceType::basis_matrix_type>::value);
  // check that the trial space reconstructs a fom state that is compatible
  // with the fom state used in the fom system class
  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value);

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

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

  // check that the trial subspace meets the right concept
  static_assert(PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value);
  // check that the fom system meets the right concept
  static_assert(FullyDiscreteSystemWithJacobianAction<
		FomSystemType,
		TotalNumberOfDesiredStates,
		typename TrialSubspaceType::basis_matrix_type>::value);
  // check that the trial space reconstructs a fom state that is compatible
  // with the fom state used in the fom system class
  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value);

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

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
#endif  // PRESSIOROM_ROM_GALERKIN_UNSTEADY_IMPLICIT_HPP_
