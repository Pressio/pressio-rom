/*
//@HEADER
// ************************************************************************
//
// qr_fwd.hpp
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

#ifndef PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_FWD_HPP_
#define PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_FWD_HPP_

namespace pressio{ namespace qr{

template<typename derived, typename matrix_t> class QRInPlaceBase;
template<typename derived, typename matrix_t, typename Q_type> class QROutOfPlaceBase;
template<typename derived, typename R_type> class RFactorBase;
template<typename derived> class QRSolveBase;

namespace impl{

#if defined PRESSIO_ENABLE_TPL_TRILINOS
#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class matrix_t, class R_t> class TpetraMVHouseholderUsingEigen;
#endif // PRESSIO_ENABLE_TPL_EIGEN

template<class matrix_t, class R_t> class TpetraMVTSQR;
template<class matrix_t, class R_t> class ModGramSchmidtMVTpetra;
template<class matrix_t, class R_t> class TpetraBlockMVTSQR;
#endif //PRESSIO_ENABLE_TPL_TRILINOS

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<typename matrix_type, typename R_t = void, typename enable=void>
class QRHouseholderDenseEigenMatrix;
#endif

template<class matrix_type, class algorithm, bool in_place, class R_type, class enable = void>
class QRSolver;

}//end namespace pressio::qr::impl


/* Note to devs: if you add a method here, and code it
 * remember to update the corresponding static_assert
 * message in the qr_traits.hpp */

struct ModifiedGramSchmidt{};
struct Householder{};

#if defined PRESSIO_ENABLE_TPL_TRILINOS
struct TSQR{};
#endif


template<class matrix_type, class algorithm, bool in_place = false>
using QRSolver = impl::QRSolver<matrix_type, algorithm, in_place, void>;

// template<class matrix_type, class algorithm, class R_type, bool in_place = false>
// using QRSolverWrapR = impl::QRSolver<matrix_type, algorithm, in_place, R_type>;

}}//end namespace pressio::qr
#endif  // PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_FWD_HPP_
