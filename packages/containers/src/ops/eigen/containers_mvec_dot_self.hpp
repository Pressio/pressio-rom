/*
//@HEADER
// ************************************************************************
//
// containers_mvec_dot_self.hpp
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
#ifndef CONTAINERS_SRC_OPS_EIGEN_MULTI_VECTOR_DOT_SELF_HPP_
#define CONTAINERS_SRC_OPS_EIGEN_MULTI_VECTOR_DOT_SELF_HPP_

namespace pressio{ namespace containers{ namespace ops{

/*
 * multi_vector dot self
 */

/* void specialize for:
 * result_t = dense dynamic eigen matrix wrapper
 */
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::details::traits<result_t>::is_dynamic
    > * = nullptr
  >
void dot_self(const mvec_t & A, result_t & C)
{
  static_assert(containers::meta::are_scalar_compatible<mvec_t, result_t>::value,
		"Types are not scalar compatible");

  const auto nAcols = A.data()->cols();
  // since C is dynamic, I can resize if needed
  if(C.data()->rows() != nAcols || C.data()->cols() != nAcols)
    C.data()->resize( nAcols, nAcols );

  *C.data() = A.data()->transpose() * (*A.data());
}


/* non void specialize for:
 * result_t = dense dynamic eigen matrix wrapper
 */
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::details::traits<result_t>::is_dynamic
    > * = nullptr
  >
result_t dot_self(const mvec_t & A)
{
  static_assert(containers::meta::are_scalar_compatible<mvec_t, result_t>::value,
		"Types are not scalar compatible");

  const auto numVecsA = A.numVectors();
  result_t C(numVecsA, numVecsA);
  dot_self(A, C);
  return C;
}


}}}//end namespace pressio::containers::ops
#endif
