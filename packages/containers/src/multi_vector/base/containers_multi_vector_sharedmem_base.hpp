/*
//@HEADER
// ************************************************************************
//
// containers_multi_vector_sharedmem_base.hpp
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

#ifndef CONTAINERS_MULTIVECTOR_BASE_MULTIVECTOR_SHAREDMEM_BASE_HPP_
#define CONTAINERS_MULTIVECTOR_BASE_MULTIVECTOR_SHAREDMEM_BASE_HPP_

#include "../containers_multi_vector_traits.hpp"

namespace pressio{ namespace containers{

template<typename derived_type>
class MultiVectorSharedMemBase
  : private utils::details::CrtpBase<
  MultiVectorSharedMemBase<derived_type>>
{
  static_assert(details::traits<derived_type>::is_shared_mem==1,
  "OOPS: distributed concrete vector inheriting from sharedMem base!");

private:
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using ord_t = typename details::traits<derived_type>::ordinal_t;

public:
  ord_t numVectors() const{
    return this->underlying().numVectorsImpl();
  }

  ord_t length() const {
    return this->underlying().lengthImpl();
  };

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_type>::type;

  using this_t = MultiVectorSharedMemBase<derived_type>;
  friend utils::details::CrtpBase<this_t>;

  MultiVectorSharedMemBase() = default;
  ~MultiVectorSharedMemBase() = default;

};//end class

}}//end namespace pressio::containers
#endif
