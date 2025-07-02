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

#ifndef PRESSIOROM_ROM_CONCEPTS_HPP_
#define PRESSIOROM_ROM_CONCEPTS_HPP_

#include "./concepts_helpers.hpp"
#include "./predicates.hpp"

namespace pressio{ namespace rom{

// -----------------------------------------------------------------------------
// Trait: ReducedState<T>
//
// Determines whether a type T is a valid reduced state type in the Pressio ROM
// framework.
//
// Specialization:
//   - When PRESSIO_ENABLE_TPL_EIGEN is defined, this trait is specialized
//     for Eigen vector types via ::pressio::is_vector_eigen<T>, and inherits
//     from std::true_type.
//
// Usage:
//   - Used in static_asserts and SFINAE to constrain templates to valid reduced
//     state vector types.
//
// Example:
//   static_assert(ReducedState<MyVectorType>::value, "Invalid reduced state type");
//
// Notes:
//   - Additional specializations can be added for other TPLs like Kokkos
//     if and when we add support for that.
// -----------------------------------------------------------------------------

template<class T, class = void>
struct ReducedState : std::false_type{};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct ReducedState<
  T, std::enable_if_t< ::pressio::is_vector_eigen<T>::value > > : std::true_type{};
#endif




// -----------------------------------------------------------------------------
// Trait: PossiblyAffineTrialColumnSubspace<T>
//
// This trait determines whether a type T conforms to the interface expected
// of a trial column subspace used in pressio ROMs.

template<class T, class enable = void>
struct PossiblyAffineTrialColumnSubspace : std::false_type{};

template<class T>
struct PossiblyAffineTrialColumnSubspace<
  T,
  std::enable_if_t<
    // must define T::reduced_state_type
    // and it  must be a valid reduced state
    ::pressio::has_reduced_state_typedef<T>::value &&
    ReducedState<typename T::reduced_state_type>::value &&

    // must define T::basis_matrix_type
    ::pressio::has_basis_matrix_typedef<T>::value &&

    // must define T::full_state_type
    ::pressio::has_full_state_typedef<T>::value &&

    // full state and basis matrix must be copy-constructible
    std::is_copy_constructible< typename T::full_state_type>::value &&
    std::is_copy_constructible< typename T::basis_matrix_type>::value &&

    // all must have valid Traits specializations
    all_have_traits<
      typename T::reduced_state_type,
      typename T::full_state_type,
      typename T::basis_matrix_type>::value &&

    // full state must be rank-1
    Traits<typename T::full_state_type>::rank == 1 &&

    // basis matrix must be rank-2
    Traits<typename T::basis_matrix_type>::rank == 2 &&

    // T must be move-constructible and move assignable
    std::is_move_constructible<T>::value &&
    std::is_assignable<T&, T&&>::value &&

    // must have: createReducedState() ->reduced_state_type
    has_const_create_reduced_state_return_result<T>::value &&

    // must have: createFullState() ->full_state_type
    has_const_create_full_state_return_result<T>::value &&

    // must have: createFullStateFromReducedState(reduced) ->full_state_type
    has_const_create_full_state_from_reduced_state<T>::value &&

    // must have: mapFromReducedState(reduced, full) ->void
    has_const_map_from_reduced_state_return_void<T>::value &&

    // must have: basisOfTranslatedSpace() ->const basis_matrix_type&
    std::is_same<
      decltype( std::declval<T const>().basisOfTranslatedSpace() ),
      typename T::basis_matrix_type const &
      >::value &&

    // must have: translationVector() ->const full_state_type&
    std::is_same<
      decltype(std::declval<T const>().translationVector()),
      const typename T::full_state_type &
      >::value &&

    // must have: basis() ->const basis_matrix_type&
    std::is_same<
        decltype( std::declval<T const>().basis() ),
        const typename T::basis_matrix_type &
      >::value &&

    // must have: dimension() ->integral type
    std::is_integral<
      decltype( std::declval<T const>().dimension() )
      >::value &&

    // must have: isColumnSpace() ->bool
    std::is_same<
      decltype( std::declval<T const>().isColumnSpace() ),
      bool
      >::value &&

    // must have: isRowSpace() ->bool
    std::is_same<
      decltype( std::declval<T const>().isRowSpace() ),
      bool
      >::value
    >
  > : std::true_type{};


// -----------------------------------------------------------------------------
// Trait: PossiblyAffineRealValuedTrialColumnSubspace<T>
//
// Builds on PossiblyAffineTrialColumnSubspace<T> by further requiring that
// all scalar types (underlying the reduced state, full state, and basis matrix)
// are floating-point types.

template<class T, class enable = void>
struct PossiblyAffineRealValuedTrialColumnSubspace: std::false_type{};

template<class T>
struct PossiblyAffineRealValuedTrialColumnSubspace<
  T,
  std::enable_if_t<
    // must already satisfy all affine subspace checks
    PossiblyAffineTrialColumnSubspace<T>::value &&

    // reduced state's scalar type must be floating-point
    std::is_floating_point< scalar_trait_t<typename T::reduced_state_type> >::value &&

    // full state's scalar type must be floating-point
    std::is_floating_point< scalar_trait_t<typename T::full_state_type> >::value &&

    // basis matrix's scalar type must be floating-point
    std::is_floating_point< scalar_trait_t<typename T::basis_matrix_type> >::value
  >
> : std::true_type{};




// -----------------------------------------------------------------------------
// Trait: SteadyFomWithJacobianAction<T, JacobianActionOperandType>
//
// This trait determines whether a type `T` represents a valid *steady* full-order
// model (FOM) interface with support for computing residuals and applying the
// Jacobian to an operand (Jacobian action).
//
// This trait inherits from std::true_type if all of the following are satisfied:
//   - `T` defines the required type aliases (state_type, residual_type).
//   - `T` provides factory methods to create residuals and Jacobian action results.
//   - `T` exposes a const-qualified method named:
//       void residualAndJacobianAction(
//         const state_type &,
//         residual_type &,
//         const JacobianActionOperandType &,
//         std::optional<fom_jac_action_t<T, JacobianActionOperandType> *>
//       )
//
// Otherwise, the trait inherits from std::false_type.
//
// This is typically used to statically constrain the type `T` in the ROM API
// when constructing a steady ROM problem.
// -----------------------------------------------------------------------------

template<class T, class JacobianActionOperandType, class enable = void>
struct SteadyFomWithJacobianAction : std::false_type{};

template<class T, class JacobianActionOperandType>
struct SteadyFomWithJacobianAction<
  T, JacobianActionOperandType,
  std::enable_if_t<

    // Must define: T::state_type
    ::pressio::has_state_typedef<T>::value &&

    // Must define: T::residual_type
    ::pressio::has_residual_typedef<T>::value &&

    // State and residual types must be copy-constructible
    std::is_copy_constructible<typename T::state_type>::value &&
    std::is_copy_constructible<typename T::residual_type>::value &&

    // Must provide: T::createResidual() returning residual_type
    ::pressio::rom::has_const_create_residual_method_return_result<
      T, typename T::residual_type >::value &&

    // Must provide: T::createResultOfJacobianActionOn(operand) returning non-void
    ::pressio::rom::has_const_create_result_of_jacobian_action_on<
      T, JacobianActionOperandType>::value &&

    // The resulting jacobian action type must be copy-constructible
    std::is_copy_constructible<
      impl::fom_jac_action_t<T, JacobianActionOperandType>
      >::value &&

    // Must implement:
    // void residualAndJacobianAction(state, residual, operand, optional<jac_action*>)
    std::is_void<
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




// -----------------------------------------------------------------------------
// Trait: SemiDiscreteFom<T>
//
// Determines whether a type `T` represents a valid semi-discrete FOM
// interface for time-dependent problems of the form:
//
//     d/dt y(t) = f(y, t)
//
// The model must provide time, state, and RHS (velocity) types, as well as
// a const-qualified method to compute the RHS given a state and time.
//
// If all requirements are met, the trait inherits from std::true_type.
// Otherwise, it inherits from std::false_type.
//
// This trait is used in Pressio to statically validate that a type conforms
// to the minimal interface needed for time integration of ODE-based FOMs.
// -----------------------------------------------------------------------------

template<class T, class enable = void>
struct SemiDiscreteFom : std::false_type{};

template<class T>
struct SemiDiscreteFom<
  T,
  std::enable_if_t<

    // Must define: T::time_type
    ::pressio::has_time_typedef<T>::value &&

    // Must define: T::state_type
    ::pressio::has_state_typedef<T>::value &&

    // Must define: T::rhs_type (i.e., velocity or RHS vector)
    ::pressio::has_rhs_typedef<T>::value &&

    // state_type must be copy-constructible
    std::is_copy_constructible<typename T::state_type>::value &&

    // rhs_type must be copy-constructible
    std::is_copy_constructible<typename T::rhs_type>::value &&

    // Must provide: createRhs() ->rhs_type
    ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::rhs_type
    >::value &&

    // Must implement: rhs(state, time, result) ->void
    ::pressio::rom::has_const_rhs_method_accept_state_indvar_result_return_void<
      T,
      typename T::state_type,
      typename T::time_type,
      typename T::rhs_type
    >::value
   >
  > : std::true_type{};

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




// -----------------------------------------------------------------------------
// Trait: SemiDiscreteFomWithJacobianAction<T, JacobianActionOperandType>
//
// Determines whether a type `T` represents a valid time-dependent full-order model (FOM)
// interface (i.e., a semi-discrete system) that also supports applying the Jacobian to
// an operand (Jacobian-action interface).
//
// Inherits from std::true_type if all the following conditions are satisfied:
//   - `T` satisfies the SemiDiscreteFom trait.
//   - `T` can create a result container for a Jacobian-action via
//         createResultOfJacobianActionOn(operand)
//   - The resulting Jacobian-action object is copy-constructible.
//   - `T` implements a const-qualified method of the form:
//         void applyJacobian(
//           const state_type&,
//           const operand_type&,
//           const time_type&,
//           result_type& )
//
// Otherwise, inherits from std::false_type.
// This is typically used to statically constrain the type `T` in the ROM API
// when constructing unsteady implicit ROM problems.
// -----------------------------------------------------------------------------

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




// -----------------------------------------------------------------------------
// Trait: SemiDiscreteFomWithMassMatrixAction<T, MassMatrixActionOperandType>
//
// Determines whether a type `T` represents a valid time-dependent full-order model (FOM)
// that also supports applying the *mass matrix* to a given operand, i.e., implements
// a mass matrix action interface.
//
// If all conditions are met, the trait inherits from std::true_type.
// Otherwise, it inherits from std::false_type.
//
// Required interface for T:
//   - Satisfies the `SemiDiscreteFom` trait.
//   - Provides:
//       - createResultOfMassMatrixActionOn(operand) ->returns result object
//       - applyMassMatrix(state, operand, time, result) ->void
//
// -----------------------------------------------------------------------------

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




// -----------------------------------------------------------------------------
// Trait: SemiDiscreteFomWithJacobianAndMassMatrixAction<T, OperandType>
//
// Determines whether a type `T` represents a valid time-dependent full-order model (FOM)
// that supports both:
//   - Jacobian-action interface (applyJacobian)
//   - Mass matrix action interface (applyMassMatrix)
//
// This trait inherits from std::true_type only if both of the following traits
// are satisfied for the same operand type:
//   - SemiDiscreteFomWithJacobianAction<T, OperandType>
//   - SemiDiscreteFomWithMassMatrixAction<T, OperandType>
//
// Otherwise, it inherits from std::false_type.
//
// This trait is used when both the Jacobian and mass matrix actions
// are required (e.g., implicit Galerkin with mass matrix).
// -----------------------------------------------------------------------------

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




// -----------------------------------------------------------------------------
// Trait: FullyDiscreteSystemWithJacobianAction<T, TotalNumStates, OperandType>
//
// Determines whether a type `T` conforms to the interface of a fully-discrete
// full-order model (FOM) that supports computing residuals and Jacobian actions
// for multiple time levels (e.g., backward differentiation formulas).
//
// This trait is enabled only if all of the following hold:
//   - The number of states (`TotalNumStates`) is either 2 or 3 (currently supported).
//   - `T` defines the following typedefs:
//       -time_type, state_type, discrete_residual_type
//   - `state_type` and `discrete_residual_type` are copy-constructible.
//   - `T` provides:
//       - createDiscreteTimeResidual() ->discrete_residual_type
//       - createResultOfDiscreteTimeJacobianActionOn(Operand) ->copyable result type
//   - `T` implements a method of the form:
//
//         void discreteTimeResidualAndJacobianAction(
//           StepCountType stepNumber,
//           time_type time,
//           std::array<state_type, TotalNumStates> const& states,
//           discrete_residual_type& residual,
//           OperandType const& operand,
//           ResultType& result);
//
// This trait is typically used to enforce constraints for FOMs employed in
// projection-based ROMs using fully-discrete time integration schemes.
// -----------------------------------------------------------------------------

template<class T, int TotalNumStates, class JacobianActionOperandType, class = void>
struct FullyDiscreteSystemWithJacobianAction : std::false_type{};

template<class T, int TotalNumStates, class JacobianActionOperandType>
struct FullyDiscreteSystemWithJacobianAction<
  T, TotalNumStates, JacobianActionOperandType,
  std::enable_if_t<

    // Supported values: must be 2 or 3 states for fully-discrete stepping
    (TotalNumStates == 2 || TotalNumStates == 3) &&

    // Must define time_type, state_type, discrete_residual_type
    ::pressio::has_time_typedef<T>::value &&
    ::pressio::has_state_typedef<T>::value &&
    ::pressio::has_discrete_residual_typedef<T>::value &&

    // state_type and residual_type must be copyable
    std::is_copy_constructible<typename T::state_type>::value &&
    std::is_copy_constructible<typename T::discrete_residual_type>::value &&

    // Must provide: createDiscreteTimeResidual() ->discrete_residual_type
    std::is_same<
      typename T::discrete_residual_type,
      decltype(std::declval<T const>().createDiscreteTimeResidual())
    >::value &&

    // Must provide: createResultOfDiscreteTimeJacobianActionOn(operand) ->copyable result type
    std::is_copy_constructible<
      decltype(
        std::declval<T const>().createResultOfDiscreteTimeJacobianActionOn(
          std::declval<JacobianActionOperandType const &>()
        )
      )
    >::value &&

    // Must implement: discreteTimeResidualAndJacobianAction(...) ->void
    ::pressio::rom::has_const_discrete_residual_jacobian_action_method<
      T,
      TotalNumStates,
      typename ::pressio::ode::StepCount::value_type,
      typename T::time_type,
      typename T::state_type,
      typename T::discrete_residual_type,
      JacobianActionOperandType,
      impl::fully_discrete_fom_jac_action_t<T, JacobianActionOperandType>
    >::value

    >
  > : std::true_type{};

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

}} // end namespace pressio::rom

#endif  // PRESSIOROM_ROM_CONCEPTS_HPP_
