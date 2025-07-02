/*
//@HEADER
// ************************************************************************
//
// qr_traits.hpp
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

#ifndef PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_TRAITS_HPP_
#define PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_TRAITS_HPP_

namespace pressio{

/* common to all cases */
template<typename matrix_type, typename algo, bool in_place>
struct qr_traits_shared_all{
  using matrix_t  = matrix_type;
  using algo_t	  = algo;
  using sc_t	  = typename ::pressio::Traits<matrix_type>::scalar_type;
  static constexpr bool in_place_ = in_place;
};


#ifdef PRESSIO_ENABLE_TPL_EIGEN
/*
 * specialize for:
 *  Eigen::DenseMatrix,
 *  R_type = void
 */
template<typename matrix_type, typename algo_t, bool in_place>
struct Traits<
  qr::impl::QRSolver<matrix_type, algo_t, in_place, void>,
    std::enable_if_t<
      ::pressio::is_dense_matrix_eigen<matrix_type>::value
      >
  > : qr_traits_shared_all<matrix_type, algo_t, in_place>
{

  static_assert( std::is_same<algo_t, qr::ModifiedGramSchmidt>::value or
     std::is_same<algo_t, qr::Householder>::value,
  "Currently, only ModifiedGramSchmidt and Householder are available for Eigen matrices.");

  using traits_all_t  = qr_traits_shared_all<matrix_type, algo_t, in_place>;
  using typename traits_all_t::matrix_t;
  using typename traits_all_t::sc_t;

  using impl_t = qr::impl::QRHouseholderDenseEigenMatrix<matrix_t, void>;
  using Q_type   = typename impl_t::Q_type;
  using concrete_t  = qr::impl::QRSolver<matrix_type, algo_t, in_place, void>;
  using inplace_base_t  = qr::QRInPlaceBase<concrete_t, matrix_type>;
  using outplace_base_t = qr::QROutOfPlaceBase<concrete_t, matrix_type, Q_type>;
  using base_compute_t  = typename std::conditional<in_place, inplace_base_t, outplace_base_t>::type;
  using base_solve_t  = qr::QRSolveBase<concrete_t>;
};
#endif


#if defined PRESSIO_ENABLE_TPL_TRILINOS

/* helpers to define the type of the implementation class to use */
template <class matrix_t, class algo_tag, class R_t, class enable = void>
struct impl_class_helper{};

template <class matrix_t, class R_t>
struct impl_class_helper<
  matrix_t, qr::TSQR, R_t,
  std::enable_if_t<
    ::pressio::is_multi_vector_tpetra<matrix_t>::value
    >
  >
{
  using impl_t = qr::impl::TpetraMVTSQR<matrix_t, R_t>;
};

template <class matrix_t, class R_t>
struct impl_class_helper<
  matrix_t, qr::TSQR, R_t,
  std::enable_if_t<
    ::pressio::is_multi_vector_tpetra_block<matrix_t>::value
    >
  >
{
  using impl_t = qr::impl::TpetraBlockMVTSQR<matrix_t, R_t>;
};

template <class matrix_t, class R_t>
struct impl_class_helper<
  matrix_t, qr::ModifiedGramSchmidt, R_t,
  std::enable_if_t<
    ::pressio::is_multi_vector_tpetra<matrix_t>::value
    >
  >
{
  using impl_t = qr::impl::ModGramSchmidtMVTpetra<matrix_t, R_t>;
};

template <class matrix_t, class R_t>
struct impl_class_helper<
  matrix_t, qr::Householder, R_t,
  std::enable_if_t<
    ::pressio::is_multi_vector_tpetra<matrix_t>::value
    >
  >{
  using impl_t = qr::impl::TpetraMVHouseholderUsingEigen<matrix_t, R_t>;
};

/*
 * specialize for Tpetra::MultiVector, R_type = void
 */
template<typename matrix_type, typename algo_t, bool in_place>
struct Traits<
  qr::impl::QRSolver<matrix_type, algo_t, in_place, void>,
    std::enable_if_t<::pressio::is_multi_vector_tpetra<matrix_type>::value>
  >
  : qr_traits_shared_all<matrix_type, algo_t, in_place>
{

  static_assert(
    std::is_same<algo_t, qr::ModifiedGramSchmidt>::value or
    std::is_same<algo_t, qr::Householder>::value or
    std::is_same<algo_t, qr::TSQR>::value,
    "Currently, only TSQR, ModifiedGramSchmidt and Householder are available for Tpetra dense matrices. \
    Use TSQR because it is fast and accurate. ModifiedGramSchmidt and Householder \
    are just here for testing purposes. ");

  using traits_all_t  = qr_traits_shared_all<matrix_type, algo_t, in_place>;
  using typename traits_all_t::matrix_t;
  using typename traits_all_t::sc_t;
  using node_t = typename matrix_type::node_type;

  using impl_t   = typename impl_class_helper<matrix_t, algo_t, void>::impl_t;
  using Q_type  = typename impl_t::Q_type;

  using concrete_t = qr::impl::QRSolver<matrix_type, algo_t, in_place, void>;
  using inplace_base_t  = qr::QRInPlaceBase<concrete_t, matrix_type>;
  using outplace_base_t = qr::QROutOfPlaceBase<concrete_t, matrix_type, Q_type>;
  using base_compute_t = typename std::conditional<in_place, inplace_base_t, outplace_base_t>::type;
  using base_solve_t = qr::QRSolveBase<concrete_t>;
};

/*
 * specialize for Tpetra::BlockMultiVector, R_type = void
 */
template<typename matrix_type, typename algo_t, bool in_place>
struct Traits<
  qr::impl::QRSolver<matrix_type, algo_t, in_place, void>,
  std::enable_if_t< ::pressio::is_multi_vector_tpetra_block<matrix_type>::value >
  >
  : qr_traits_shared_all<matrix_type, algo_t, in_place>
{

  static_assert( std::is_same<algo_t, qr::TSQR>::value,
    "Currently, only TSQR is available for BlockTpetra dense matrices.");

  using traits_all_t  = qr_traits_shared_all<matrix_type, algo_t, in_place>;
  using typename traits_all_t::matrix_t;
  using typename traits_all_t::sc_t;
  using node_t = typename matrix_type::node_type;

  using impl_t   = typename impl_class_helper<matrix_t, algo_t, void>::impl_t;
  using Q_type  = typename impl_t::Q_type;

  using concrete_t = qr::impl::QRSolver<matrix_type, algo_t, in_place, void>;
  using inplace_base_t  = qr::QRInPlaceBase<concrete_t, matrix_type>;
  using outplace_base_t = qr::QROutOfPlaceBase<concrete_t, matrix_type, Q_type>;
  using base_compute_t = typename std::conditional<in_place,inplace_base_t, outplace_base_t>::type;
  using base_solve_t = qr::QRSolveBase<concrete_t>;
};
#endif // PRESSIO_ENABLE_TPL_TRILINOS

}
#endif  // PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_TRAITS_HPP_
