/*
//@HEADER
// ************************************************************************
//
// pressio_containers_multi_vector_include.hpp
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

#ifndef CONTAINERS_MULTI_VECTOR_PRESSIO_CONTAINERS_MULTI_VECTOR_INCLUDE_HPP_
#define CONTAINERS_MULTI_VECTOR_PRESSIO_CONTAINERS_MULTI_VECTOR_INCLUDE_HPP_

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./predicates/containers_native_kokkos_multi_vector_meta.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./predicates/containers_native_epetra_multi_vector_meta.hpp"
#include "./predicates/containers_native_tpetra_block_multi_vector_meta.hpp"
#include "./predicates/containers_native_tpetra_multi_vector_meta.hpp"
#endif
#include "./predicates/containers_native_eigen_multi_vector_meta.hpp"
#include "./predicates/containers_native_arbitrary_multi_vector_meta.hpp"


#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./predicates/containers_is_multi_vector_wrapper_epetra.hpp"
#include "./predicates/containers_is_multi_vector_wrapper_tpetra.hpp"
#include "./predicates/containers_is_multi_vector_wrapper_tpetra_block.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./predicates/containers_is_multi_vector_wrapper_kokkos.hpp"
#endif
#include "./predicates/containers_is_multi_vector_wrapper_eigen.hpp"
#include "./predicates/containers_is_multi_vector_wrapper_arbitrary.hpp"
#include "./predicates/containers_is_multi_vector_wrapper.hpp"

#include "./containers_multi_vector_traits.hpp"


#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./concrete/containers_multi_vector_distributed_epetra.hpp"
#include "./concrete/containers_multi_vector_distributed_tpetra_block.hpp"
#include "./concrete/containers_multi_vector_distributed_tpetra.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./concrete/containers_multi_vector_sharedmem_kokkos.hpp"
#endif
#include "./concrete/containers_multi_vector_arbitrary.hpp"
#include "./concrete/containers_multi_vector_sharedmem_eigen_dynamic.hpp"



#endif  // CONTAINERS_MULTI_VECTOR_PRESSIO_CONTAINERS_MULTI_VECTOR_INCLUDE_HPP_
