/*
//@HEADER
// ************************************************************************
//
// solvers_linear_solver_selector_impl.hpp
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

#ifndef PRESSIOROM_SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_SOLVER_SELECTOR_IMPL_HPP_
#define PRESSIOROM_SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_SOLVER_SELECTOR_IMPL_HPP_

#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "solvers_linear_eigen_direct_impl.hpp"
#include "solvers_linear_eigen_iterative_impl.hpp"
#include "solvers_linear_eigen_iterative_matrix_free_impl.hpp"
#endif

// because this uses teuchos lapack wrapper
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "solvers_linear_kokkos_direct_getrs_impl.hpp"
#endif

namespace pressio{ namespace linsol{ namespace impl{


template<typename TagType, typename MatrixType, typename enable = void>
struct Selector{
  using type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<typename UserDefinedLinearOperatorType>
struct Selector<
  iterative::GMRES, UserDefinedLinearOperatorType, void
  >
{
  using tag_t = iterative::GMRES;
  using solver_traits = ::pressio::linsol::Traits<tag_t>;
  using type = EigenIterativeMatrixFreeWrapper<tag_t, UserDefinedLinearOperatorType>;
};

template<typename TagType, typename MatrixType>
struct Selector<
  TagType, MatrixType,
  std::enable_if_t<
    ::pressio::linsol::Traits<TagType>::iterative and
    (::pressio::is_dense_matrix_eigen<MatrixType>::value or
     ::pressio::is_sparse_matrix_eigen<MatrixType>::value)
    >
  >
{
  using solver_traits = ::pressio::linsol::Traits<TagType>;
  using type = ::pressio::linsol::impl::EigenIterativeWrapper<TagType, MatrixType>;
};

template<typename TagType, typename MatrixType>
struct Selector<
  TagType, MatrixType,
  std::enable_if_t<
    ::pressio::linsol::Traits<TagType>::direct and
    (::pressio::is_dense_matrix_eigen<MatrixType>::value or
     ::pressio::is_sparse_matrix_eigen<MatrixType>::value)>
  >
{
  using solver_traits = ::pressio::linsol::Traits<TagType>;
  using type = ::pressio::linsol::impl::EigenDirectWrapper<TagType, MatrixType>;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template<typename TagType, typename MatrixType>
struct Selector<
  TagType, MatrixType,
  std::enable_if_t<
    ::pressio::linsol::Traits<TagType>::direct and
    (::pressio::is_dense_matrix_kokkos<MatrixType>::value)
    >
  >
{
  using solver_traits   = ::pressio::linsol::Traits<TagType>;
  using type = ::pressio::linsol::impl::KokkosDirectGETRS<MatrixType>;
};
#endif

}}}
#endif  // PRESSIOROM_SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_SOLVER_SELECTOR_IMPL_HPP_
