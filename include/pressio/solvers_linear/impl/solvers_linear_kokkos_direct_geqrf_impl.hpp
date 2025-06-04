/*
//@HEADER
// ************************************************************************
//
// solvers_linear_kokkos_direct_geqrf_impl.hpp
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

#ifndef PRESSIO_SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_KOKKOS_DIRECT_GEQRF_IMPL_HPP_
#define PRESSIO_SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_KOKKOS_DIRECT_GEQRF_IMPL_HPP_

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialDenseSolver.hpp>
#include <Teuchos_BLAS.hpp>
#endif

#if defined PRESSIO_ENABLE_TPL_KOKKOS and defined KOKKOS_ENABLE_CUDA
#include <cuda_runtime.h>
#include <cusolverDn.h>
#endif

namespace pressio { namespace linearsolvers{ namespace impl{

template<typename SolverT, typename MatrixType, typename enable = void>
class KokkosDirect;

template<typename MatrixType>
class KokkosDirect<::pressio::linearsolvers::direct::geqrf, MatrixType>
{
public:
  using solver_tag	= ::pressio::linearsolvers::direct::geqrf;
  using this_type          = KokkosDirect<solver_tag, MatrixType>;
  using matrix_type	= MatrixType;
  using scalar_type        = typename MatrixType::value_type;
  using solver_traits   = ::pressio::linearsolvers::Traits<solver_tag>;

  static_assert( solver_traits::kokkos_enabled == true,
  		 "the native solver must suppport kokkos to use in KokkosDirect");
  static_assert( solver_traits::direct == true,
  		 "the native solver must be direct to use in KokkosDirect");

public:
  KokkosDirect(){
#ifdef KOKKOS_ENABLE_CUDA
    auto cusolverStatus = cusolverDnCreate(&cuDnHandle_);
    assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
#endif
  }

  KokkosDirect(const KokkosDirect &) = delete;

  ~KokkosDirect(){
#ifdef KOKKOS_ENABLE_CUDA
    auto cusolverStatus = cusolverDnDestroy(cuDnHandle_);
    assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
#endif
  }

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  /*
   * enable if:
   * the matrix has layout left (i.e. column major)
   * T is a kokkos vector
   * T has host execution space
   * T and MatrixType have same execution space
   */
  template <typename _MatrixType = MatrixType, typename T>
  std::enable_if_t<
    std::is_same<typename _MatrixType::traits::array_layout, Kokkos::LayoutLeft>::value
    and ::pressio::is_vector_kokkos<T>::value
    /*and ::pressio::containers::details::traits<T>::has_host_execution_space and*/
    and std::is_same<typename T::traits::execution_space, typename _MatrixType::traits::execution_space>::value
  >
  solve(const _MatrixType & A, const T& b, T & y)
  {
    const auto Aext0 = A.extent(0);
    const auto Aext1 = A.extent(1);
    if (Aext0 != auxMat_.extent(0) or	Aext1 != auxMat_.extent(1)){
    	Kokkos::resize(auxMat_, Aext0, Aext1);
    }

    Kokkos::deep_copy(auxMat_, A);
    this->solveAllowMatOverwrite(auxMat_, b, y);
  }

private:
  /*
   * enable if:
   * the matrix has layout left (i.e. column major)
   * T is a kokkos vector
   * T has host execution space
   * T and MatrixType have same execution space
   */
  template <typename _MatrixType = MatrixType, typename T>
  std::enable_if_t<
    std::is_same<typename _MatrixType::traits::array_layout, Kokkos::LayoutLeft>::value
    and ::pressio::is_vector_kokkos<T>::value
    /*::pressio::containers::details::traits<T>::has_host_execution_space and*/
    and std::is_same<typename T::traits::execution_space, typename _MatrixType::traits::execution_space>::value
  >
  solveAllowMatOverwrite(_MatrixType & A, const T& b, T & y)
  {
    // gerts is for square matrices
    assert(A.extent(0) == A.extent(1));

    assert(A.extent(0) == b.extent(0) );
    assert(A.extent(1) == y.extent(0) );

    // only one rhs
    constexpr int nRhs = 1;
    // just use n, since rows == cols
    const auto n = A.extent(0);

    // to store the return code of the function
    int info = 0;

    if (tau_.size() != n){
      tau_.resize(n);
    }

    if (lwork_ == -1){
      // this means we need to query what lwork should be
      // and resize the work array "work_" properly
      lpk_.GEQRF(n, n, A.data(), n, tau_.data(), work_.data(), lwork_, &info);
      work_.resize(work_[0]);
      lwork_ = work_[0];
    }
    // do QR
    lpk_.GEQRF(n, n, A.data(), n, tau_.data(), work_.data(), lwork_, &info);
    assert(info == 0);
    // std::cout << " info " << info
    // 	      << " lwork " << lwork_
    // 	      << " workSz " << work_.size() << std::endl;

    // Note that A is overwritten starting from here!
    // we need to deep copy b into y and pass y to ormqr
    // because it is overwritten with Q^T b
    Kokkos::deep_copy(y, b);

    // compute Q^T b
    const char side = 'L';
    const char trans = 'T';
    lpk_.ORMQR(side, trans, n, nRhs, n,
    	       A.data(),
    	       n,
    	       tau_.data(),
    	       y.data(),
    	       n,
    	       work_.data(),
    	       lwork_,
    	       &info);
    //std::cout << " info-ormqr " << info << std::endl;

    // solver R y = Q^T b
    constexpr scalar_type alpha = static_cast<scalar_type>(1);
    blas_.TRSM(Teuchos::ESide::LEFT_SIDE,
    	       Teuchos::EUplo::UPPER_TRI,
    	       Teuchos::ETransp::NO_TRANS,
    	       Teuchos::EDiag::NON_UNIT_DIAG,
    	       n, nRhs, alpha,
    	       A.data(), n, y.data(), n);
  }
#endif

private:
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  Teuchos::LAPACK<int, scalar_type> lpk_;
  Teuchos::BLAS<int, scalar_type> blas_;

  // lwork is used for the geqrf, but if you read the doc,
  // if lwork == -1, then geqrf does a query of what is needed.
  // more details are shown in the code above on how we use lwork
  int lwork_ = -1;
  std::vector<scalar_type> work_ = {0};
  std::vector<scalar_type> tau_ = {};
  MatrixType auxMat_ = {};
#endif

#if defined PRESSIO_ENABLE_TPL_KOKKOS and defined KOKKOS_ENABLE_CUDA
  cusolverDnHandle_t cuDnHandle_;
#endif
};

}}} // end namespace pressio::linearsolvers::impl
#endif  // PRESSIO_SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_KOKKOS_DIRECT_GEQRF_IMPL_HPP_
