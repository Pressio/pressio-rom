/*
//@HEADER
// ************************************************************************
//
// qr_in_place_base.hpp
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

#ifndef PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_BASE_CLASSES_HPP_
#define PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_BASE_CLASSES_HPP_

namespace pressio{ namespace qr{

template <typename T, typename Q_T, typename Enable = void>
struct is_legitimate_vector_type_for_qr_project : std::false_type {};

template <typename T, typename Q_t>
struct is_legitimate_vector_type_for_qr_project<T, Q_t,
   std::enable_if_t<
       ::pressio::Traits<T>::rank == 1 and
     // the vector type should be from same package as Q
     (false
#ifdef PRESSIO_ENABLE_TPL_EIGEN
     or (::pressio::is_native_container_eigen<T>::value and ::pressio::is_native_container_eigen<Q_t>::value)
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
     or (::pressio::is_native_container_kokkos<T>::value and ::pressio::is_native_container_kokkos<Q_t>::value)
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
     or (::pressio::is_vector_tpetra<T>::value and ::pressio::is_multi_vector_tpetra<Q_t>::value)
     or (::pressio::is_vector_tpetra_block<T>::value and ::pressio::is_multi_vector_tpetra_block<Q_t>::value)
#endif
     )
   >
> : std::true_type{};

template<typename DerivedType, typename MatrixType>
class QRInPlaceBase
{

  using this_t = QRInPlaceBase<DerivedType, MatrixType>;

  /* workaround for nvcc issue with templates,
  see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<DerivedType>::type;

public:
  void computeThin(MatrixType & A){
    static_cast<DerivedType &>(*this).computeThinImpl(A);
  }

private:
  QRInPlaceBase() = default;
  ~QRInPlaceBase() = default;
};

template<typename DerivedType, typename MatrixType, typename QType>
class QROutOfPlaceBase
{

  using this_t = QROutOfPlaceBase<DerivedType, MatrixType, QType>;

  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<DerivedType>::type;

public:
  void computeThin(const MatrixType & A){
    static_cast<DerivedType &>(*this).computeThinImpl(A);
  }

  const QType & cRefQFactor() const {
    return static_cast<DerivedType const &>(*this).cRefQFactorImpl();
  }

  template <typename VectorInType, typename VectorOutType>
  std::enable_if_t<
    ::pressio::Traits<VectorInType>::rank ==1 and
    ::pressio::Traits<VectorOutType>::rank ==1 and
    is_legitimate_vector_type_for_qr_project<VectorInType, QType>::value
  >
  applyQTranspose(const VectorInType & vecIn, VectorOutType & vecOut) const{
    static_cast<DerivedType const &>(*this).applyQTransposeImpl(vecIn, vecOut);
  }

  template <typename VectorInType, typename VectorOutType>
  std::enable_if_t<
    ::pressio::Traits<VectorInType>::rank ==1 and
    ::pressio::Traits<VectorOutType>::rank ==1
  >
  applyRTranspose(const VectorInType & vecIn, VectorOutType & vecOut) const{
    static_cast<DerivedType const &>(*this).applyRTransposeImpl(vecIn, vecOut);
  }

  QROutOfPlaceBase() = default;
  ~QROutOfPlaceBase() = default;
};


template<typename DerivedType, typename RType>
class RFactorBase
{

  using this_t = RFactorBase<DerivedType, RType>;

  /* workaround for nvcc issue with templates, 
  see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<DerivedType>::type;

public:
  const RType & cRefRFactor() const {
    return static_cast<DerivedType const&>(*this).cRefRFactorImpl();
  }

  RFactorBase() = default;
  ~RFactorBase() = default;

  mutable std::shared_ptr<RType> Rmat_ = nullptr;
};


template<typename DerivedType>
class QRSolveBase
{
  using this_t = QRSolveBase<DerivedType>;

  /* workaround for nvcc issue with templates, 
  see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<DerivedType>::type;

public:
  template <typename VectorType>
  std::enable_if_t<
   ::pressio::Traits<VectorType>::rank ==1 
  >
  solve(const VectorType & rhs, VectorType & y)const {
    static_cast<DerivedType const &>(*this).solveImpl(rhs, y);
  }

  QRSolveBase() = default;
  ~QRSolveBase() = default;
};

}}//end namespace pressio::qr
#endif  // PRESSIOROM_SOLVERS_NONLINEAR_IMPL_QR_QR_BASE_CLASSES_HPP_
