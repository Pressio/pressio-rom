/*
//@HEADER
// ************************************************************************
//
// galerkin_unsteady_system_default_rhs_with_mass_matrix.hpp
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

#ifndef PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_DEFAULT_RHS_WITH_MASS_MATRIX_HPP_
#define PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_DEFAULT_RHS_WITH_MASS_MATRIX_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class ReducedMassMatType,
  class TrialSubspaceType,
  class FomSystemType
  >
class GalerkinDefaultOdeSystemOnlyRhsAndMassMatrix
{
  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;
  using fom_mm_action_result_type =
    decltype(std::declval<FomSystemType const>().createResultOfMassMatrixActionOn
	     (std::declval<basis_matrix_type const &>()) );

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using rhs_type                  = ReducedRhsType;
  using mass_matrix_type          = ReducedMassMatType;

  GalerkinDefaultOdeSystemOnlyRhsAndMassMatrix(const TrialSubspaceType & trialSubspace,
					       const FomSystemType & fomSystem)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      fomRhs_(fomSystem.createRhs()),
      fomMMAction_(fomSystem.createResultOfMassMatrixActionOn(trialSubspace_.get().basisOfTranslatedSpace()))
  {}

  GalerkinDefaultOdeSystemOnlyRhsAndMassMatrix(GalerkinDefaultOdeSystemOnlyRhsAndMassMatrix const &) = delete;
  GalerkinDefaultOdeSystemOnlyRhsAndMassMatrix& operator=(GalerkinDefaultOdeSystemOnlyRhsAndMassMatrix const&) = delete;
  GalerkinDefaultOdeSystemOnlyRhsAndMassMatrix(GalerkinDefaultOdeSystemOnlyRhsAndMassMatrix &&) = default;
  GalerkinDefaultOdeSystemOnlyRhsAndMassMatrix& operator=(GalerkinDefaultOdeSystemOnlyRhsAndMassMatrix &&) = default;

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  rhs_type createRhs() const{
    return impl::CreateGalerkinRhs<rhs_type>()(trialSubspace_.get().dimension());
  }

  mass_matrix_type createMassMatrix() const{
    return impl::CreateGalerkinMassMatrix<mass_matrix_type>()(trialSubspace_.get().dimension());
  }

  void massMatrixAndRhs(const state_type & reducedState,
			const IndVarType & rhsEvaluationTime,
			mass_matrix_type & reducedMassMat,
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

    fomSystem_.get().applyMassMatrix(fomState_, phi, rhsEvaluationTime, fomMMAction_);

    // compute the reduced mass matrix
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
			    alpha, phi, fomMMAction_,
			    beta, reducedMassMat);
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  mutable typename FomSystemType::rhs_type fomRhs_;
  mutable fom_mm_action_result_type fomMMAction_;
};

}}} // end pressio::rom::impl
#endif  // PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_DEFAULT_RHS_WITH_MASS_MATRIX_HPP_
