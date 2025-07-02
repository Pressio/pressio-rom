/*
//@HEADER
// ************************************************************************
//
// lspg_steady_system_default.hpp
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

#ifndef PRESSIOROM_ROM_IMPL_LSPG_STEADY_SYSTEM_DEFAULT_HPP_
#define PRESSIOROM_ROM_IMPL_LSPG_STEADY_SYSTEM_DEFAULT_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
LSPG steady default represents:

  min_x ||fom_r(phi x)||

- fom_r is the fom "residual"
- phi is the basis
*/
template <
  class ReducedStateType,
  class TrialSubspaceType,
  class FomSystemType,
  class PossiblyRefWrapperOperatorScalerType
  >
class LspgSteadyDefaultSystem
{

  // need to deduce the type of the action of the fom jacobian
  // which becomes the jacobian_type of the problem
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createResultOfJacobianActionOn
	     (std::declval<typename TrialSubspaceType::basis_matrix_type const &>())
	     );

public:
  // aliases required by the pressio solvers
  using state_type    = ReducedStateType;
  using residual_type = typename FomSystemType::residual_type;
  using jacobian_type = fom_jac_action_result_type;

  // here _RawScalerType must be a template because it is the raw scaler type
  // which can be forwarded to the reference wrapper potentially
  template<class _RawScalerType>
  LspgSteadyDefaultSystem(const TrialSubspaceType & trialSubspace,
			  const FomSystemType & fomSystem,
			  _RawScalerType && scaler)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      scaler_(std::forward<_RawScalerType>(scaler))
  {}

  LspgSteadyDefaultSystem(LspgSteadyDefaultSystem const &) = delete;
  LspgSteadyDefaultSystem& operator=(LspgSteadyDefaultSystem const&) = delete;
  LspgSteadyDefaultSystem(LspgSteadyDefaultSystem &&) = default;
  LspgSteadyDefaultSystem& operator=(LspgSteadyDefaultSystem &&) = default;

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  residual_type createResidual() const{
    return fomSystem_.get().createResidual();
  }

  jacobian_type createJacobian() const{
    return fomSystem_.get().createResultOfJacobianActionOn
      (trialSubspace_.get().basisOfTranslatedSpace());
  }

  void residualAndJacobian(const state_type & lspgState,
			   residual_type & lspgResidual,
			   std::optional<jacobian_type *> lspgJacobian) const
  {
    trialSubspace_.get().mapFromReducedState(lspgState, fomState_);

    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    fomSystem_.get().residualAndJacobianAction(fomState_,
					       lspgResidual,
					       phi, lspgJacobian);
    scaler_(fomState_, lspgResidual, lspgJacobian);
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  PossiblyRefWrapperOperatorScalerType scaler_;
};


}}} // end pressio::rom::impl
#endif  // PRESSIOROM_ROM_IMPL_LSPG_STEADY_SYSTEM_DEFAULT_HPP_
