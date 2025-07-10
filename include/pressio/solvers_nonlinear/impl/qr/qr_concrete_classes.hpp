/*
//@HEADER
// ************************************************************************
//
// qr_concrete_classes.hpp
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

#ifndef PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_CONCRETE_CLASSES_HPP_
#define PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_CONCRETE_CLASSES_HPP_

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include <Teuchos_SerialQRDenseSolver.hpp>
#endif


namespace pressio{ namespace qr{ namespace impl{

/* specialize for R_type == void */
template<typename MatrixType, typename AlgoType>
class QRSolver<MatrixType, AlgoType, true, void> 
  : public ::pressio::Traits< QRSolver<MatrixType, AlgoType,true, void>>::base_compute_t,
    public ::pressio::Traits< QRSolver<MatrixType, AlgoType,true, void>>::base_solve_t
{

  static_assert(::pressio::Traits<MatrixType>::rank == 2, 
    "QRSolver only supports rank-2 objects");

  using this_t		= QRSolver<MatrixType, AlgoType, true, void>;
  using qr_traits	= ::pressio::Traits<this_t>;
  using base_compute_t	= typename qr_traits::base_compute_t;
  using base_solve_t	= typename qr_traits::base_solve_t;

  using impl_t	       = typename qr_traits::impl_t;
  impl_t myImpl_;

  void computeThinImpl(MatrixType & A){
    myImpl_.computeThinInPlace(A);
  }

  template <typename VectorType>
  void solveImpl(const VectorType & rhs, VectorType & y) const {
    myImpl_.template doLinSolve<VectorType>(rhs, y);
  }

public:
  QRSolver() = default;
  ~QRSolver() = default;

private:
  friend base_compute_t;
  friend base_solve_t;
};


/* R_type == void, in_place = false */
template<typename MatrixType, typename algo>
class QRSolver<MatrixType, algo, false, void>
  : public ::pressio::Traits< QRSolver<MatrixType, algo, false, void>>::base_compute_t,
    public ::pressio::Traits< QRSolver<MatrixType, algo, false, void>>::base_solve_t
{

  static_assert(::pressio::Traits<MatrixType>::rank == 2,
    "QRSolver only supports rank-2 objects");

  using this_t         = QRSolver<MatrixType, algo, false, void>;
  using qr_traits      = ::pressio::Traits<this_t>;
  using base_compute_t = typename qr_traits::base_compute_t;
  using base_solve_t   = typename qr_traits::base_solve_t;
  using Q_t            = typename qr_traits::Q_type;
  using impl_t         = typename qr_traits::impl_t;
  impl_t myImpl_;

public:
  QRSolver() = default;
  ~QRSolver() = default;

private:
  void computeThinImpl(const MatrixType & A){
    myImpl_.computeThinOutOfPlace(A);
  }

  template < typename VectorInType, typename VectorOutType>
  void applyQTransposeImpl(const VectorInType & vecIn, VectorOutType & vecOut) const{
    myImpl_.template applyQTranspose<VectorInType, VectorOutType>(vecIn, vecOut);
  }

  template < typename VectorInType, typename VectorOutType>
  void applyRTransposeImpl(const VectorInType & vecIn, VectorOutType & vecOut) const{
    myImpl_.template applyRTranspose<VectorInType, VectorOutType>(vecIn, vecOut);
  }

  template <typename VectorType>
  void solveImpl(const VectorType & rhs, VectorType & y)const{
    myImpl_.template doLinSolve<VectorType>(rhs, y);
  }

  const Q_t & cRefQFactorImpl() const {
    return myImpl_.QFactor();
  }

private:
  friend base_compute_t;
  friend base_solve_t;
};


}}} // end namespace pressio::qr::impl
#endif  // PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_CONCRETE_CLASSES_HPP_
