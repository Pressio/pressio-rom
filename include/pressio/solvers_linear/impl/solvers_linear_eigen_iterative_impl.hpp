/*
//@HEADER
// ************************************************************************
//
// solvers_linear_eigen_iterative_impl.hpp
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

#ifndef PRESSIO_SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_EIGEN_ITERATIVE_IMPL_HPP_
#define PRESSIO_SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_EIGEN_ITERATIVE_IMPL_HPP_

#include "solvers_linear_iterative_base.hpp"

namespace pressio { namespace linearsolvers{ namespace impl{

template<typename TagType, typename MatrixType>
class EigenIterativeWrapper
  : public IterativeBase< EigenIterativeWrapper<TagType, MatrixType>>
{
  using this_type = EigenIterativeWrapper<TagType, MatrixType>;
  using solver_traits = ::pressio::linearsolvers::Traits<TagType>;
  using native_solver_type = typename solver_traits::template eigen_solver_type<MatrixType>;
  using base_iterative_type = IterativeBase<this_type>;
  using iteration_type = typename base_iterative_type::iteration_type;

  static_assert( solver_traits::eigen_enabled == true,
		 "the native solver must be from Eigen to use in EigenIterativeWrapper");
  static_assert( solver_traits::direct == false,
		 "The native eigen solver must be iterative to use in EigenIterativeWrapper");

public:
  using matrix_type    = MatrixType;
  using scalar_type    = typename MatrixType::Scalar;

  EigenIterativeWrapper() = default;
  EigenIterativeWrapper(EigenIterativeWrapper const &) = delete;
  EigenIterativeWrapper& operator=(EigenIterativeWrapper const &) = delete;

  iteration_type numIterationsExecuted() const{
    return mysolver_.iterations();
  }

  scalar_type finalError() const{
    return mysolver_.error();
  }

  template <typename T>
  void solve(const MatrixType & A, const T& b, T & y){
    this->resetLinearSystem(A);
    mysolver_.setMaxIterations(this->maxIters_);
    y = mysolver_.solve(b);
  }

private:
  void resetLinearSystem(const MatrixType& A){
    mysolver_.setMaxIterations(this->maxIters_);
    mysolver_.compute(A);
  }

  friend base_iterative_type;
  native_solver_type mysolver_ = {};
};

}}} // end namespace pressio::solvers::iterarive::impl
#endif  // PRESSIO_SOLVERS_LINEAR_IMPL_SOLVERS_LINEAR_EIGEN_ITERATIVE_IMPL_HPP_
