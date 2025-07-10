/*
//@HEADER
// ************************************************************************
//
// reduced_operators_traits.hpp
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

#ifndef PRESSIOROM_ROM_REDUCED_OPERATORS_TRAITS_HPP_
#define PRESSIOROM_ROM_REDUCED_OPERATORS_TRAITS_HPP_

namespace pressio{ namespace rom{

/*
  steady galerkin
*/
template<class T, class = void>
struct SteadyGalerkinDefaultReducedOperatorsTraits
{
  using reduced_state_type    = void;
  using reduced_residual_type = void;
  using reduced_jacobian_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct SteadyGalerkinDefaultReducedOperatorsTraits<
  T, std::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using reduced_state_type    = T;
  using reduced_residual_type = T;

  // figure out what is the reduced jacobian type
  // if the reduced state is Eigen vector,
  // it makes sense to use an Eigen dense matrix to store
  // the Galerkin jacobian since all reduced operators are dense
  using reduced_jacobian_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
};
#endif

/*
  unsteady explicit galerkin
*/
template<class T, class = void>
struct ExplicitGalerkinDefaultReducedOperatorsTraits
{
  using reduced_state_type    = void;
  using reduced_rhs_type = void;
  using reduced_mass_matrix_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct ExplicitGalerkinDefaultReducedOperatorsTraits<
  T, std::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using reduced_state_type    = T;
  using reduced_rhs_type = T;
  using reduced_mass_matrix_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
};
#endif


/*
  unsteady implicit galerkin
*/
template<class T, class = void>
struct ImplicitGalerkinDefaultReducedOperatorsTraits
{
  using reduced_state_type    = void;
  using reduced_residual_type = void;
  using reduced_jacobian_type = void;
  using reduced_mass_matrix_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct ImplicitGalerkinDefaultReducedOperatorsTraits<
  T, std::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using reduced_state_type    = T;
  using reduced_residual_type = T;
  using reduced_jacobian_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
  using reduced_mass_matrix_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
};
#endif


/*
  steady LSPG
*/
template<class T, class = void>
struct SteadyLspgDefaultReducedOperatorsTraits
{
  using reduced_state_type = void;
  using hessian_type    = void;
  using gradient_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct SteadyLspgDefaultReducedOperatorsTraits<
  T, std::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using reduced_state_type = T;
  using gradient_type = T;
  using hessian_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
};
#endif

/*
  unsteady LSPG
*/
template<class T, class = void>
struct UnsteadyLspgDefaultReducedOperatorsTraits
{
  using reduced_state_type = void;
  using hessian_type    = void;
  using gradient_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct UnsteadyLspgDefaultReducedOperatorsTraits<
  T, std::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using reduced_state_type = T;
  using gradient_type = T;
  using hessian_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
};
#endif

}}
#endif  // PRESSIOROM_ROM_REDUCED_OPERATORS_TRAITS_HPP_
