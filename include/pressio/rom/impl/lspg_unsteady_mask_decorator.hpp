/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_cont_time_decorators.hpp
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

#ifndef PRESSIO_ROM_IMPL_LSPG_UNSTEADY_MASK_DECORATOR_HPP_
#define PRESSIO_ROM_IMPL_LSPG_UNSTEADY_MASK_DECORATOR_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class MaskerType,
  class ResidualType,
  class JacobianType,
  class Maskable
  >
class LspgMaskDecorator : public Maskable
{
  using ind_var_type = typename Maskable::independent_variable_type;
  using unmasked_residual_type = typename Maskable::residual_type;
  using unmasked_jacobian_type = typename Maskable::jacobian_type;

public:
  // required
  using independent_variable_type = typename Maskable::independent_variable_type;
  using state_type    = typename Maskable::state_type;
  using residual_type = ResidualType;
  using jacobian_type = JacobianType;

  LspgMaskDecorator() = delete;

  template <class TrialSubspaceType, class FomSystemType>
  LspgMaskDecorator(const TrialSubspaceType & trialSubspace,
		    const FomSystemType & fomSystem,
		    FomStatesManager<TrialSubspaceType> & fomStatesManager,
		    const MaskerType & masker)
    : Maskable(trialSubspace, fomSystem, fomStatesManager),
      masker_(masker),
      unMaskedResidual_(Maskable::createResidual()),
      unMaskedJacobian_(Maskable::createJacobian())
  {}

public:
  state_type createState() const{
    return Maskable::createState();
  }

  residual_type createResidual() const{
    auto tmp = Maskable::createResidual();
    return masker_.get().createResultOfMaskActionOn(tmp);
  }

  jacobian_type createJacobian() const{
    auto tmp = Maskable::createJacobian();
    return masker_.get().createResultOfMaskActionOn(tmp);
  }

  template <class StencilStatesContainerType, class StencilRhsContainerType>
  void operator()(::pressio::ode::StepScheme odeSchemeName,
		  const state_type & predictedReducedState,
		  const StencilStatesContainerType & reducedStatesStencilManager,
		  StencilRhsContainerType & fomRhsStencilManger,
		  const ::pressio::ode::StepEndAt<ind_var_type> & rhsEvaluationTime,
		  ::pressio::ode::StepCount step,
		  const ::pressio::ode::StepSize<ind_var_type> & dt,
		  residual_type & R,
		  std::optional<jacobian_type *> Jo) const
  {
    std::optional<unmasked_jacobian_type*> unmaskedJo{&unMaskedJacobian_};

    Maskable::operator()(odeSchemeName, predictedReducedState,
			 reducedStatesStencilManager, fomRhsStencilManger,
			 rhsEvaluationTime, step, dt,
			 unMaskedResidual_, unmaskedJo);
    masker_(unMaskedResidual_, R);

    if (Jo){
     masker_(unMaskedJacobian_, *Jo.value());
    }
  }

private:
  std::reference_wrapper<const MaskerType> masker_;
  mutable unmasked_residual_type unMaskedResidual_;
  mutable unmasked_jacobian_type unMaskedJacobian_;
};

}}}
#endif  // PRESSIO_ROM_IMPL_LSPG_UNSTEADY_MASK_DECORATOR_HPP_
