/*
//@HEADER
// ************************************************************************
//
// lspg_helpers.hpp
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

#ifndef PRESSIOROM_ROM_IMPL_LSPG_HELPERS_HPP_
#define PRESSIOROM_ROM_IMPL_LSPG_HELPERS_HPP_

namespace pressio{ namespace rom{ namespace impl{

template<class T= void>
void valid_scheme_for_lspg_else_throw(::pressio::ode::StepScheme name){
  if (   name != ::pressio::ode::StepScheme::BDF1
      && name != ::pressio::ode::StepScheme::BDF2)
  {
    throw std::runtime_error("LSPG currently accepting BDF1 or BDF2");
  }
}

template<class TrialSubspaceType, class FomSystemType>
void lspg_static_check_trial_and_system(const TrialSubspaceType & trialSpace,
					const FomSystemType & fomSystem)
{
  static_assert(PossiblyAffineTrialColumnSubspace<TrialSubspaceType>::value,
     "You are trying to create a steady lspg problem but the \
trialSpace does not meet the required PossiblyAffineTrialColumnSubspace concept.");

  static_assert
    (SteadyFomWithJacobianAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
     "You are trying to create a steady lspg problem but the \
FOM system does not meet the required SteadyFomWithJacobianAction concept.");
}

template<typename T>
void steady_lspg_static_check_api_return_type(){
  static constexpr bool val =
#ifdef PRESSIO_ENABLE_CXX20
    nonlinearsolvers::NonlinearSystemFusingResidualAndJacobian<T>;
#else
  nonlinearsolvers::NonlinearSystemFusingResidualAndJacobian<T>::value;
#endif
  static_assert(val,
		"The return type must satisify the NonlinearSystemFusingResidualAndJacobian concept.");
}

}}} // end pressio::rom::impl
#endif  // PRESSIOROM_ROM_IMPL_LSPG_HELPERS_HPP_
