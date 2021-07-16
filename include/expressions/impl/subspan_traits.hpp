/*
//@HEADER
// ************************************************************************
//
// containers_subspan_traits.hpp
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

#ifndef CONTAINERS_EXPRESSIONS_SUBSPAN_CONTAINERS_SUBSPAN_TRAITS_HPP_
#define CONTAINERS_EXPRESSIONS_SUBSPAN_CONTAINERS_SUBSPAN_TRAITS_HPP_

namespace pressio{ namespace expressions{ namespace impl{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename MatrixType>
struct subspan_traits<
  SubspanExpr<MatrixType>,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dense_matrix_eigen<
     typename std::remove_cv<MatrixType>::type
    >::value
    >
  >
  : public containers_shared_traits<PackageIdentifier::Eigen, true, 2>,
  public matrix_shared_traits<traits<MatrixType>::is_sparse>
{
  static constexpr bool is_static = true;
  static constexpr bool is_dynamic  = !is_static;

  using scalar_type  = typename traits<MatrixType>::scalar_type;
  using ordinal_type = typename traits<MatrixType>::ordinal_type;
  using size_type    = ordinal_type;

  // the reference type is conditionnal because the native expression
  // returns by value when object is const
  using reference_type = typename std::conditional<
    std::is_const<MatrixType>::value, scalar_type, scalar_type &
  >::type;

  using const_reference_type = typename std::conditional<
    std::is_const<MatrixType>::value, scalar_type, scalar_type const &
  >::type;

  // type of the native expression
  using _native_expr_type = decltype(
    std::declval<MatrixType>().block( std::declval<size_type>(),
				     std::declval<size_type>(),
                                     std::declval<size_type>(),
				     std::declval<size_type>() )
    );
  using _const_native_expr_type = decltype(
    std::declval<const MatrixType>().block( std::declval<size_type>(),
					   std::declval<size_type>(),
                                           std::declval<size_type>(),
					   std::declval<size_type>() )
    );
  using native_expr_type = typename std::conditional<
    std::is_const<MatrixType>::value,
    _const_native_expr_type,
    _native_expr_type
  >::type;

  using const_data_return_type = native_expr_type const *;
  using data_return_type = native_expr_type *;
};
#endif


// #ifdef PRESSIO_ENABLE_TPL_KOKKOS
// template <typename matrix_type>
// struct traits<
//   ::pressio::expressions::SubspanExpr<matrix_type>,
//   ::pressio::mpl::enable_if_t<
//     ::pressio::is_dense_matrix_kokkos<
//      typename std::remove_cv<matrix_type>::type
//     >::value
//     >
//   >
//   : public containers_shared_traits<PackageIdentifier::Kokkos, true, 2>,
//   public matrix_shared_traits<false>
// {
//   static constexpr bool is_static = true;
//   static constexpr bool is_dynamic  = !is_static;

//   // using wrapped_t = typename traits<matrix_type>::wrapped_t;
//   using scalar_t	= typename traits<matrix_type>::scalar_t;
//   using execution_space = typename traits<matrix_type>::execution_space;
//   using memory_space	= typename traits<matrix_type>::memory_space;
//   using device_t	= typename traits<matrix_type>::device_t;
//   using device_type	= typename traits<matrix_type>::device_t;
//   using ordinal_t	= typename traits<matrix_type>::ordinal_t;
//   using size_t		= ordinal_t;
//   using pair_t		= std::pair<size_t, size_t>;
//   using reference_t	  = scalar_t &;
//   using const_reference_t = scalar_t const &;

//   using _native_expr_t = decltype
//     (
//      Kokkos::subview(std::declval<matrix_type>(),
// 		     std::declval<pair_t>(),
// 		     std::declval<pair_t>())
//     );
//   using _const_native_expr_t = decltype
//     (
//      Kokkos::subview(std::declval<const matrix_type>(),
// 		     std::declval<pair_t>(),
// 		     std::declval<pair_t>())
//      );
//   using native_expr_t = typename std::conditional<
//     std::is_const<matrix_type>::value,
//     _const_native_expr_t,
//     _native_expr_t
//   >::type;

//   using const_data_return_t = native_expr_t const *;
//   using data_return_t	    = native_expr_t *;
// };
// #endif

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template <typename matrix_type>
// struct traits<
//   ::pressio::containers::expressions::SubspanExpr<matrix_type>,
//   ::pressio::mpl::enable_if_t<
//     ::pressio::containers::predicates::is_rank2_tensor_wrapper_pybind<matrix_type>::value
//     >
//   >
//   : public containers_shared_traits<
//   typename details::traits<matrix_type>::wrapped_t,
//   WrappedPackageIdentifier::Pybind, true, 2>,
//   public matrix_shared_traits<false>
// {
//   static constexpr bool is_static = true;
//   static constexpr bool is_dynamic  = !is_static;

//   using wrapped_t = typename traits<matrix_type>::wrapped_t;
//   using scalar_t  = typename traits<matrix_type>::scalar_t;
//   using ordinal_t = typename traits<matrix_type>::ordinal_t;
//   using size_t    = ordinal_t;
//   using reference_t =  scalar_t &;
//   using const_reference_t = scalar_t const &;
// };
// #endif

}}}
#endif  // CONTAINERS_EXPRESSIONS_SUBSPAN_CONTAINERS_SUBSPAN_TRAITS_HPP_
