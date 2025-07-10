/*
//@HEADER
// ************************************************************************
//
// galerkin_unsteady_system_masked_rhs_only.hpp
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

#ifndef PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_ONLY_HPP_
#define PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_ONLY_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  masked explicit galerkin system represents:

     d hat{y}/dt = HrOp mask( fom_rhs(phi*hat{y}, ...) )

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- HrOp is the hyper-red operator
*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType
  >
class GalerkinMaskedOdeSystemOnlyRhs
{
  // deduce types
  using unmasked_fom_rhs_type = typename FomSystemType::rhs_type;
  using masked_fom_rhs_type =
    decltype(std::declval<MaskerType const>().createResultOfMaskActionOn
	     (std::declval<unmasked_fom_rhs_type const &>()));

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using rhs_type      = ReducedRhsType;

  GalerkinMaskedOdeSystemOnlyRhs(const TrialSubspaceType & trialSubspace,
				 const FomSystemType & fomSystem,
				 const MaskerType & masker,
				 const HyperReducerType & hyperReducer)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      hyperReducer_(hyperReducer),
      masker_(masker),
      unMaskedFomRhs_(fomSystem.createRhs()),
      maskedFomRhs_(masker.createResultOfMaskActionOn(unMaskedFomRhs_))
  {}


  GalerkinMaskedOdeSystemOnlyRhs(GalerkinMaskedOdeSystemOnlyRhs const &) = delete;
  GalerkinMaskedOdeSystemOnlyRhs& operator=(GalerkinMaskedOdeSystemOnlyRhs const&) = delete;
  GalerkinMaskedOdeSystemOnlyRhs(GalerkinMaskedOdeSystemOnlyRhs &&) = default;
  GalerkinMaskedOdeSystemOnlyRhs& operator=(GalerkinMaskedOdeSystemOnlyRhs &&) = default;

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  rhs_type createRhs() const{
    return impl::CreateGalerkinRhs<rhs_type>()(trialSubspace_.get().dimension());
  }

  void rhs(const state_type & reducedState,
	   const IndVarType & rhsEvaluationTime,
	   rhs_type & reducedRhs) const
  {

    // reconstruct fom state fomState = phi*reducedState
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);
    // evaluate fomRhs
    fomSystem_.get().rhs(fomState_, rhsEvaluationTime, unMaskedFomRhs_);
    // mask it
    masker_(unMaskedFomRhs_, maskedFomRhs_);
    // evaluate reduced rhs
    hyperReducer_(maskedFomRhs_, rhsEvaluationTime, reducedRhs);
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReducerType> hyperReducer_;
  std::reference_wrapper<const MaskerType> masker_;
  mutable unmasked_fom_rhs_type unMaskedFomRhs_;
  mutable masked_fom_rhs_type maskedFomRhs_;
};

}}} // end pressio::rom::impl
#endif  // PRESSIOROM_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_ONLY_HPP_
