/*
//@HEADER
// ************************************************************************
//
// qr.hpp
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

#ifndef PRESSIO_SOLVERS_NONLINEAR_IMPL_QR_HPP_
#define PRESSIO_SOLVERS_NONLINEAR_IMPL_QR_HPP_

#include "qr/qr_fwd.hpp"
#include "qr/qr_base_classes.hpp"
#include "qr/qr_traits.hpp"
#include "qr/qr_concrete_classes.hpp"

#ifdef PRESSIO_ENABLE_TPL_EIGEN
#include "qr/qr_eigen_impl.hpp"
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "qr/qr_tpetra_impl.hpp"
#ifdef PRESSIO_ENABLE_EPETRA
#include "qr/qr_epetra_multi_vector_tsqr_impl.hpp"
#include "qr/qr_epetra_mv_householder_using_eigen_impl.hpp"
#include "qr/qr_epetra_multi_vector_modified_gram_schmidt_impl.hpp"
#endif // PRESSIO_ENABLE_EPETRA
#endif // PRESSIO_ENABLE_TPL_TRILINOS

#endif  // PRESSIO_SOLVERS_NONLINEAR_IMPL_QR_HPP_
