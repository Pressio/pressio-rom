/*
//@HEADER
// ************************************************************************
//
// rom_custom_ops_for_linear_decoder.hpp
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

#ifndef ROM_WILL_BE_CONCEPTS_CUSTOM_OPS_ROM_CUSTOM_OPS_FOR_LINEAR_DECODER_HPP_
#define ROM_WILL_BE_CONCEPTS_CUSTOM_OPS_ROM_CUSTOM_OPS_FOR_LINEAR_DECODER_HPP_

namespace pressio{ namespace rom{ namespace concepts {

template<
  typename T,
  typename mat_type, typename rom_state_type, typename fom_state_type,
  typename enable = void
  >
struct custom_ops_for_linear_decoder
  : std::false_type{};


template <
  typename T, typename mat_type, typename rom_state_type, typename fom_state_type
  >
struct custom_ops_for_linear_decoder<
  T, mat_type, rom_state_type, fom_state_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::ops::predicates::has_void_method_product_mat_vec<
      T,
      ::pressio::nontranspose,
      typename ::pressio::containers::details::traits<mat_type>::scalar_t,
      typename ::pressio::containers::details::traits<mat_type>::wrapped_t,
      rom_state_type,
      typename ::pressio::containers::details::traits<fom_state_type>::wrapped_t
      >::value and
    ::pressio::ops::predicates::has_void_method_product_mat_vec<
      T,
      ::pressio::nontranspose,
      typename ::pressio::containers::details::traits<mat_type>::scalar_t,
      typename ::pressio::containers::details::traits<mat_type>::wrapped_t,
      typename ::pressio::containers::details::traits<rom_state_type>::span_const_ret_t,
      typename ::pressio::containers::details::traits<fom_state_type>::wrapped_t
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::rom::concepts
#endif  // ROM_WILL_BE_CONCEPTS_CUSTOM_OPS_ROM_CUSTOM_OPS_FOR_LINEAR_DECODER_HPP_
