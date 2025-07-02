/*
//@HEADER
// ************************************************************************
//
// qr_eigen_dense_out_of_place_impl.hpp
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

#ifndef PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_EIGEN_IMPL_HPP_
#define PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_EIGEN_IMPL_HPP_

#include <Eigen/QR>

namespace pressio{ namespace qr{ namespace impl{

template< typename MatrixType, typename R_t>
class QRHouseholderDenseEigenMatrix<
  MatrixType, R_t,
  std::enable_if_t<
    ::pressio::is_dense_matrix_eigen<MatrixType>::value
    >
  >
{
public:
  using sc_t	       = typename ::pressio::Traits<MatrixType>::scalar_type;
  using Q_type = Eigen::Matrix<sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  using factorizer_t = Eigen::HouseholderQR<MatrixType>;

private:
  mutable std::shared_ptr<Q_type> Qmat_	     = {};
  mutable std::shared_ptr<factorizer_t> fct_ = {};

public:
  QRHouseholderDenseEigenMatrix() = default;
  ~QRHouseholderDenseEigenMatrix() = default;

  void computeThinOutOfPlace(const MatrixType & A)
  {
    const auto rows = A.rows();
    const auto cols = A.cols();
    fct_ = std::make_shared<factorizer_t>(A);

    if (!Qmat_ or (Qmat_->rows()!=rows and Qmat_->cols()!=cols)){
      Qmat_ = std::make_shared<Q_type>(rows,cols);
      Qmat_->setZero();
    }

    *Qmat_ = fct_->householderQ() * Q_type::Identity(rows,cols);
  }

  template < typename VectorInType, typename VectorOutType>
  void applyQTranspose(const VectorInType & vecIn, VectorOutType & vecOut) const
  {
    constexpr auto beta  = static_cast<sc_t>(0);
    constexpr auto alpha = static_cast<sc_t>(1);
    ::pressio::ops::product(::pressio::transpose(), alpha, *this->Qmat_, vecIn, beta, vecOut);
  }

  template < typename VectorInType, typename VectorOutType>
  void applyRTranspose(const VectorInType & vecIn, VectorOutType & y) const
  {
    // y = R^T vecIn
    auto vecSize = ::pressio::ops::extent(y, 0);
    auto & Rm = fct_->matrixQR().block(0,0,vecSize,vecSize).template triangularView<Eigen::Upper>();
    y = Rm.transpose() * vecIn;
  }

  template <typename VectorType>
  void doLinSolve(const VectorType & rhs, VectorType & y)const
  {
    auto vecSize = ::pressio::ops::extent(y, 0);
    auto & Rm = fct_->matrixQR().block(0,0,vecSize,vecSize).
      template triangularView<Eigen::Upper>();
    y = Rm.solve(rhs);
  }

  const Q_type & QFactor() const {
    return *this->Qmat_;
  }
};

}}} // end namespace pressio::qr::impl
#endif  // PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_EIGEN_IMPL_HPP_
