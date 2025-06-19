/*
//@HEADER
// ************************************************************************
//
// rom_fom_system_continuous_time.hpp
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

#ifndef PRESSIO_ROM_ROM_CONCEPTS_HPP_
#define PRESSIO_ROM_ROM_CONCEPTS_HPP_

#include "concepts_helpers.hpp"
#include "predicates.hpp"

namespace pressio{ namespace rom{

template<class T, class = void>
struct ReducedState : std::false_type{};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct ReducedState<
  T, std::enable_if_t< ::pressio::is_vector_eigen<T>::value > > : std::true_type{};
#endif

// ----------------------------------------------------------------------------
// SUBSPACES
// ----------------------------------------------------------------------------
template<class T, class enable = void>
struct PossiblyAffineTrialColumnSubspace : std::false_type{};

template<class T>
struct PossiblyAffineTrialColumnSubspace<
  T,
  std::enable_if_t<
    ::pressio::has_reduced_state_typedef<T>::value
    && ReducedState<typename T::reduced_state_type>::value
    && ::pressio::has_basis_matrix_typedef<T>::value
    && ::pressio::has_full_state_typedef<T>::value
    && std::is_copy_constructible< typename T::full_state_type>::value
    && std::is_copy_constructible< typename T::basis_matrix_type>::value
    && all_have_traits<
      typename T::reduced_state_type,
      typename T::full_state_type,
      typename T::basis_matrix_type>::value
    && Traits<typename T::full_state_type>::rank == 1
    && Traits<typename T::basis_matrix_type>::rank == 2
    /**/
    && std::is_move_constructible<T>::value
    && std::is_assignable<T&, T&&>::value
    //
    && has_const_create_reduced_state_return_result<T>::value
    && has_const_create_full_state_return_result<T>::value
    && has_const_create_full_state_from_reduced_state<T>::value
    /**/
    && has_const_map_from_reduced_state_return_void<T>::value
    && std::is_same<
      decltype( std::declval<T const>().basisOfTranslatedSpace() ),
      typename T::basis_matrix_type const &
      >::value
    && std::is_same<
      decltype(std::declval<T const>().translationVector()),
      const typename T::full_state_type &
      >::value
    && std::is_same<
        decltype( std::declval<T const>().basis() ),
        const typename T::basis_matrix_type &
      >::value
    && std::is_integral<
      decltype( std::declval<T const>().dimension() )
      >::value
    //
    && std::is_same<
      decltype( std::declval<T const>().isColumnSpace() ),
      bool
      >::value
    && std::is_same<
      decltype( std::declval<T const>().isRowSpace() ),
      bool
      >::value
    >
  > : std::true_type{};


template<class T, class enable = void>
struct PossiblyAffineRealValuedTrialColumnSubspace: std::false_type{};

template<class T>
struct PossiblyAffineRealValuedTrialColumnSubspace<
  T,
  std::enable_if_t<
  PossiblyAffineTrialColumnSubspace<T>::value
  && std::is_floating_point< scalar_trait_t<typename T::reduced_state_type> >::value
  && std::is_floating_point< scalar_trait_t<typename T::full_state_type> >::value
  && std::is_floating_point< scalar_trait_t<typename T::basis_matrix_type> >::value
  >
> : std::true_type{};


// ----------------------------------------------------------------------------
// VARIOUS
// ----------------------------------------------------------------------------
template <class T, class MaskerType, class = void>
struct MaskableWith : std::false_type{};

template <class T, class MaskerType>
struct MaskableWith<
  T, MaskerType,
  std::enable_if_t<
    std::is_copy_constructible<
      decltype
      (std::declval<MaskerType const>().createResultOfMaskActionOn
       (std::declval<T const &>())
       )
    >::value
    && std::is_void<
	decltype
	(
	 std::declval<MaskerType const>()
	 (
	   std::declval<T const &>(),
	   std::declval<impl::mask_action_t<MaskerType, T> &>()
	 )
	)
       >::value
  >
> : std::true_type{};

// ----------------------------------------------------------------------------
// SYSTEMS
// ----------------------------------------------------------------------------

template<class T, class JacobianActionOperandType, class enable = void>
struct SteadyFomWithJacobianAction : std::false_type{};

template<class T, class JacobianActionOperandType>
struct SteadyFomWithJacobianAction<
  T, JacobianActionOperandType,
  std::enable_if_t<
       ::pressio::has_state_typedef<T>::value
    && ::pressio::has_residual_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::residual_type>::value
    /**/
    && ::pressio::rom::has_const_create_residual_method_return_result<
	 T, typename T::residual_type >::value
    && ::pressio::rom::has_const_create_result_of_jacobian_action_on<
	 T, JacobianActionOperandType>::value
    && std::is_copy_constructible<
	 impl::fom_jac_action_t<T, JacobianActionOperandType>
	 >::value
    /**/
    && std::is_void<
      decltype
      (
       std::declval<T const>().residualAndJacobianAction
       (
	std::declval<typename T::state_type const&>(),
	std::declval<typename T::residual_type &>(),
	std::declval<JacobianActionOperandType const&>(),
	std::declval< std::optional<impl::fom_jac_action_t<T, JacobianActionOperandType> *> >()
	)
       )
    >::value
   >
  > : std::true_type{};

template<class T, class enable = void>
struct SemiDiscreteFom : std::false_type{};

template<class T>
struct SemiDiscreteFom<
  T,
  std::enable_if_t<
       ::pressio::has_time_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_rhs_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::rhs_type>::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::rhs_type >::value
    && ::pressio::rom::has_const_rhs_method_accept_state_indvar_result_return_void<
      T, typename T::state_type, typename T::time_type, typename T::rhs_type>::value
   >
  > : std::true_type{};


template<class T, class MassMatrixActionOperandType, class enable = void>
struct SemiDiscreteFomWithMassMatrixAction : std::false_type{};

template<class T, class MassMatrixActionOperandType>
struct SemiDiscreteFomWithMassMatrixAction<
  T, MassMatrixActionOperandType,
  std::enable_if_t<
    SemiDiscreteFom<T>::value
    //
    && std::is_copy_constructible<
      decltype
      (
       std::declval<T const>().createResultOfMassMatrixActionOn
       (
	std::declval<MassMatrixActionOperandType const &>()
	)
       )
      >::value
    && std::is_void<
       decltype
       (
	std::declval<T const>().applyMassMatrix
	(
	 std::declval<typename T::state_type const&>(),
	 std::declval<MassMatrixActionOperandType const&>(),
	 std::declval<typename T::time_type const &>(),
	 std::declval<impl::fom_mass_matrix_action_t<T,  MassMatrixActionOperandType> &>()
	 )
	)
       >::value
   >
  > : std::true_type{};


template<class T, class JacobianActionOperandType, class enable = void>
struct SemiDiscreteFomWithJacobianAction : std::false_type{};

template<class T,  class JacobianActionOperandType>
struct SemiDiscreteFomWithJacobianAction<
  T,  JacobianActionOperandType,
  std::enable_if_t<
       SemiDiscreteFom<T>::value
    && ::pressio::rom::has_const_create_result_of_jacobian_action_on<
	 T,  JacobianActionOperandType>::value
    && std::is_copy_constructible<
	 impl::fom_jac_action_t<T, JacobianActionOperandType>
	 >::value
    && ::pressio::rom::has_const_apply_jacobian_method_accept_state_operand_time_result_return_void<
	 T, typename T::state_type,
         JacobianActionOperandType, typename T::time_type,
	 impl::fom_jac_action_t<T, JacobianActionOperandType>
	 >::value
    >
  > : std::true_type{};


template<class T, class OperandType, class enable = void>
struct SemiDiscreteFomWithJacobianAndMassMatrixAction : std::false_type{};

template<class T, class OperandType>
struct SemiDiscreteFomWithJacobianAndMassMatrixAction<
  T,  OperandType,
  std::enable_if_t<
     SemiDiscreteFomWithJacobianAction<T, OperandType>::value
  && SemiDiscreteFomWithMassMatrixAction<T, OperandType>::value
  >
  > : std::true_type{};


template<class T, int TotalNumStates, class JacobianActionOperandType, class = void>
struct FullyDiscreteSystemWithJacobianAction : std::false_type{};

template<class T, int TotalNumStates, class JacobianActionOperandType>
struct FullyDiscreteSystemWithJacobianAction<
  T, TotalNumStates, JacobianActionOperandType,
  std::enable_if_t<
    (TotalNumStates == 2 || TotalNumStates == 3)
    && ::pressio::has_time_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_discrete_residual_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::discrete_residual_type>::value
    && std::is_same<
	 typename T::discrete_residual_type,
	 decltype(std::declval<T const>().createDiscreteTimeResidual())
	 >::value
    && std::is_copy_constructible<
	decltype
	(
	 std::declval<T const>().createResultOfDiscreteTimeJacobianActionOn
	  (
	   std::declval<JacobianActionOperandType const &>()
	  )
	 )
	>::value
    && ::pressio::rom::has_const_discrete_residual_jacobian_action_method<
      T, TotalNumStates,
      typename ::pressio::ode::StepCount::value_type,
      typename T::time_type,
      typename T::state_type,
      typename T::discrete_residual_type,
      JacobianActionOperandType,
      impl::fully_discrete_fom_jac_action_t<T, JacobianActionOperandType>
      >::value
    >
  > : std::true_type{};


// ----------------------------------------------------------------------------
// SYSTEMS REFINEMENTS FOR REAL VALUED
// ----------------------------------------------------------------------------

template<class T, class JacobianActionOperandType, class enable = void>
struct RealValuedSteadyFomWithJacobianAction : std::false_type{};

template<class T, class JacobianActionOperandType>
struct RealValuedSteadyFomWithJacobianAction<
  T, JacobianActionOperandType,
  std::enable_if_t<
       SteadyFomWithJacobianAction<T, JacobianActionOperandType>::value
    && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
    && std::is_floating_point< scalar_trait_t<typename T::residual_type> >::value
    && std::is_floating_point<
       scalar_trait_t< impl::fom_jac_action_t<T, JacobianActionOperandType> >
     >::value
    >
  > : std::true_type{};

// --------------------------------------------------------
template<class T, class enable = void>
struct RealValuedSemiDiscreteFom : std::false_type{};

template<class T>
struct RealValuedSemiDiscreteFom<
  T,
  std::enable_if_t<
       SemiDiscreteFom<T>::value
    && std::is_floating_point< typename T::time_type>::value
    && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
    && std::is_floating_point< scalar_trait_t<typename T::rhs_type> >::value
    >
  > : std::true_type{};

// --------------------------------------------------------
template<class T, class MassMatrixActionOperandType, class enable = void>
struct RealValuedSemiDiscreteFomWithMassMatrixAction : std::false_type{};

template<class T, class MassMatrixActionOperandType>
struct RealValuedSemiDiscreteFomWithMassMatrixAction<
  T, MassMatrixActionOperandType,
  std::enable_if_t<
       RealValuedSemiDiscreteFom<T>::value
    && SemiDiscreteFomWithMassMatrixAction<T, MassMatrixActionOperandType>::value
    && std::is_floating_point<
	 scalar_trait_t< impl::fom_mass_matrix_action_t<T, MassMatrixActionOperandType> >
      >::value
    >
  > : std::true_type{};

// --------------------------------------------------------
template<class T, class OperandType, class enable = void>
struct RealValuedSemiDiscreteFomWithJacobianAction : std::false_type{};

template<class T, class OperandType>
struct RealValuedSemiDiscreteFomWithJacobianAction<
  T, OperandType,
  std::enable_if_t<
       RealValuedSemiDiscreteFom<T>::value
    && SemiDiscreteFomWithJacobianAction<T, OperandType>::value
    && std::is_floating_point<
	 scalar_trait_t< impl::fom_jac_action_t<T, OperandType> >
        >::value
    >
  > : std::true_type{};

// --------------------------------------------------------
template<class T, class OperandType, class enable = void>
struct RealValuedSemiDiscreteFomWithJacobianAndMassMatrixAction : std::false_type{};

template<class T, class OperandType>
struct RealValuedSemiDiscreteFomWithJacobianAndMassMatrixAction<
  T, OperandType,
  std::enable_if_t<
       RealValuedSemiDiscreteFomWithJacobianAction<T, OperandType>::value
    && RealValuedSemiDiscreteFomWithMassMatrixAction<T, OperandType>::value
    >
  > : std::true_type{};

// --------------------------------------------------------
template<class T, int TotalNumStates, class OperandType, class enable = void>
struct RealValuedFullyDiscreteSystemWithJacobianAction : std::false_type{};

template<class T, int TotalNumStates, class OperandType>
struct RealValuedFullyDiscreteSystemWithJacobianAction<
  T, TotalNumStates, OperandType,
  std::enable_if_t<
    FullyDiscreteSystemWithJacobianAction<T, TotalNumStates, OperandType>::value
    && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
    && std::is_floating_point< scalar_trait_t<typename T::discrete_residual_type> >::value
    && std::is_floating_point<
      scalar_trait_t< impl::fully_discrete_fom_jac_action_t<T, OperandType> >
     >::value
    >
  > : std::true_type{};

}} // end namespace pressio::rom

#endif  // PRESSIO_ROM_ROM_CONCEPTS_HPP_
