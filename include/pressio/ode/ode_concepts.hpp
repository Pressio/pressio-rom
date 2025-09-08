/*
//@HEADER
// ************************************************************************
//
// ode_concepts.hpp
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

#ifndef PRESSIOROM_ODE_ODE_CONCEPTS_HPP_
#define PRESSIOROM_ODE_ODE_CONCEPTS_HPP_

#include "./ode_predicates.hpp"

namespace pressio{ namespace ode{

#ifdef PRESSIO_ENABLE_CXX20
template <class T>
concept OdeSystem =
  requires(){
  typename T::independent_variable_type; }
  && std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::rhs_type>
  && requires(const T & A,
	      const typename T::state_type & state,
	      const typename T::independent_variable_type & evalValue,
	      typename T::rhs_type & f)
  {
    { A.createState() } -> std::same_as<typename T::state_type>;
    { A.createRhs()   } -> std::same_as<typename T::rhs_type>;
    { A.rhs(state, evalValue, f) } -> std::same_as<void>;
  };

template <class T>
concept RealValuedOdeSystem =
     OdeSystem<T>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::rhs_type> >
  && std::convertible_to<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >;

#else

template<class T, class enable = void>
struct OdeSystem : std::false_type{};

template<class T>
struct OdeSystem<
  T,
  std::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_rhs_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::rhs_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::rhs_type >::value
    && std::is_void<
      decltype(
	       std::declval<T const>().rhs
	       (
		std::declval<typename T::state_type const&>(),
		std::declval<typename T::independent_variable_type const &>(),
		std::declval<typename T::rhs_type &>()
	       )
	   )
      >::value
   >
  > : std::true_type{};

template<class T, class enable = void>
struct RealValuedOdeSystem : std::false_type{};

template<class T>
struct RealValuedOdeSystem<
  T, std::enable_if_t<
       OdeSystem<T>::value
       && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
       && std::is_floating_point< scalar_trait_t<typename T::rhs_type> >::value
       && std::is_convertible<
	 typename T::independent_variable_type,
	 scalar_trait_t<typename T::state_type> >::value
  >
  > : std::true_type{};

#endif


#ifdef PRESSIO_ENABLE_CXX20
template <class T>
concept OdeSystemFusingRhsAndJacobian =
  requires(){ typename T::independent_variable_type; }
  && std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::rhs_type>
  && std::copy_constructible<typename T::jacobian_type>
  && requires(const T & A,
	      const typename T::state_type & state,
	      const typename T::independent_variable_type & evalValue,
	      typename T::rhs_type & f,
	      std::optional<typename T::jacobian_type*> J)
  {
    { A.createState()    } -> std::same_as<typename T::state_type>;
    { A.createRhs()      } -> std::same_as<typename T::rhs_type>;
    { A.createJacobian() } -> std::same_as<typename T::jacobian_type>;
    { A.rhsAndJacobian(state, evalValue, f, J) } -> std::same_as<void>;
  };

template <class T>
concept RealValuedOdeSystemFusingRhsAndJacobian =
     OdeSystemFusingRhsAndJacobian<T>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::rhs_type> >
  && std::floating_point< scalar_trait_t<typename T::jacobian_type> >
  && std::convertible_to<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >;

#else

template<class T, class enable = void>
struct OdeSystemFusingRhsAndJacobian : std::false_type{};

template<class T>
struct OdeSystemFusingRhsAndJacobian<
  T,
  std::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_rhs_typedef<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::rhs_type>::value
    && std::is_copy_constructible<typename T::jacobian_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::rhs_type >::value
    && ::pressio::ode::has_const_create_jacobian_method_return_result<
      T, typename T::jacobian_type >::value
    && std::is_void<
      decltype(
	       std::declval<T const>().rhsAndJacobian
	       (
		std::declval<typename T::state_type const&>(),
		std::declval<typename T::independent_variable_type const &>(),
                std::declval<typename T::rhs_type &>(),
		std::declval< std::optional<typename T::jacobian_type*> >()
	       )
	   )
      >::value
   >
  > : std::true_type{};

template<class T, class enable = void>
struct RealValuedOdeSystemFusingRhsAndJacobian : std::false_type{};

template<class T>
struct RealValuedOdeSystemFusingRhsAndJacobian<
  T,
  std::enable_if_t<
    OdeSystemFusingRhsAndJacobian<T>::value
  && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
  && std::is_floating_point< scalar_trait_t<typename T::rhs_type> >::value
  && std::is_floating_point< scalar_trait_t<typename T::jacobian_type> >::value
  && std::is_convertible<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >::value
  > > : std::true_type{};

#endif

#ifdef PRESSIO_ENABLE_CXX20
template <class T>
concept OdeSystemFusingMassMatrixAndRhs =
  requires(){ typename T::independent_variable_type; }
  && std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::rhs_type>
  && std::copy_constructible<typename T::mass_matrix_type>
  && requires(const T & A,
	      const typename T::state_type & state,
	      const typename T::independent_variable_type & evalValue,
	      typename T::mass_matrix_type & M,
	      typename T::rhs_type & f)
  {
    { A.createState()      } -> std::same_as<typename T::state_type>;
    { A.createRhs()        } -> std::same_as<typename T::rhs_type>;
    { A.createMassMatrix() } -> std::same_as<typename T::mass_matrix_type>;
    { A.massMatrixAndRhs(state, evalValue, M, f) } -> std::same_as<void>;
  };

template <class T>
concept RealValuedOdeSystemFusingMassMatrixAndRhs =
  OdeSystemFusingMassMatrixAndRhs<T>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::rhs_type> >
  && std::floating_point< scalar_trait_t<typename T::mass_matrix_type> >
  && std::convertible_to<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >;

#else

template<class T, class enable = void>
struct OdeSystemFusingMassMatrixAndRhs : std::false_type{};

template<class T>
struct OdeSystemFusingMassMatrixAndRhs<
  T,
  std::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_rhs_typedef<T>::value
    && ::pressio::has_mass_matrix_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::rhs_type>::value
    && std::is_copy_constructible<typename T::mass_matrix_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::rhs_type >::value
    && ::pressio::ode::has_const_create_mass_matrix_method_return_result<
      T, typename T::mass_matrix_type >::value
    && std::is_void<
      decltype(
	       std::declval<T const>().massMatrixAndRhs
	       (
		std::declval<typename T::state_type const&>(),
		std::declval<typename T::independent_variable_type const &>(),
		std::declval<typename T::mass_matrix_type &>(),
                std::declval<typename T::rhs_type &>()
	       )
	   )
      >::value
   >
  > : std::true_type{};

template<class T, class enable = void>
struct RealValuedOdeSystemFusingMassMatrixAndRhs : std::false_type{};

template<class T>
struct RealValuedOdeSystemFusingMassMatrixAndRhs<
  T, std::enable_if_t<
    OdeSystemFusingMassMatrixAndRhs<T>::value
  && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
  && std::is_floating_point< scalar_trait_t<typename T::rhs_type> >::value
  && std::is_floating_point< scalar_trait_t<typename T::mass_matrix_type> >::value
  && std::is_convertible<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >::value
  >
  > : std::true_type{};

#endif


#ifdef PRESSIO_ENABLE_CXX20
template <class T>
concept CompleteOdeSystem =
  requires(){ typename T::independent_variable_type; }
  && std::copy_constructible<typename T::state_type>
  && std::copy_constructible<typename T::rhs_type>
  && std::copy_constructible<typename T::mass_matrix_type>
  && std::copy_constructible<typename T::jacobian_type>
  && requires(const T & A,
	      const typename T::state_type & state,
	      const typename T::independent_variable_type & evalValue,
	      typename T::mass_matrix_type & M,
	      typename T::rhs_type & f,
	      std::optional<typename T::jacobian_type*> J)
  {
    { A.createState()      } -> std::same_as<typename T::state_type>;
    { A.createRhs()        } -> std::same_as<typename T::rhs_type>;
    { A.createMassMatrix() } -> std::same_as<typename T::mass_matrix_type>;
    { A.massMatrixAndRhsAndJacobian(state, evalValue, M, f, J) } -> std::same_as<void>;
  };

template <class T>
concept RealValuedCompleteOdeSystem =
     CompleteOdeSystem<T>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::rhs_type> >
  && std::floating_point< scalar_trait_t<typename T::mass_matrix_type> >
  && std::floating_point< scalar_trait_t<typename T::jacobian_type> >
  && std::convertible_to<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >;

#else

template<class T, class enable = void>
struct CompleteOdeSystem : std::false_type{};

template<class T>
struct CompleteOdeSystem<
  T,
  std::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_rhs_typedef<T>::value
    && ::pressio::has_mass_matrix_typedef<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::rhs_type>::value
    && std::is_copy_constructible<typename T::mass_matrix_type>::value
    && std::is_copy_constructible<typename T::jacobian_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_rhs_method_return_result<
      T, typename T::rhs_type >::value
    && ::pressio::ode::has_const_create_mass_matrix_method_return_result<
      T, typename T::mass_matrix_type >::value
    && ::pressio::ode::has_const_create_jacobian_method_return_result<
      T, typename T::jacobian_type >::value
    //
    && std::is_void<
      decltype(
	       std::declval<T const>().massMatrixAndRhsAndJacobian
	       (
		std::declval<typename T::state_type const&>(),
		std::declval<typename T::independent_variable_type const &>(),
                std::declval<typename T::mass_matrix_type &>(),
                std::declval<typename T::rhs_type &>(),
		std::declval< std::optional<typename T::jacobian_type*> >()
	       )
	   )
      >::value
   >
  > : std::true_type{};

template<class T, class enable = void>
struct RealValuedCompleteOdeSystem : std::false_type{};

template<class T>
struct RealValuedCompleteOdeSystem<
  T, std::enable_if_t<
       CompleteOdeSystem<T>::value
       && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
       && std::is_floating_point< scalar_trait_t<typename T::rhs_type> >::value
       && std::is_floating_point< scalar_trait_t<typename T::mass_matrix_type> >::value
       && std::is_floating_point< scalar_trait_t<typename T::jacobian_type> >::value
       && std::is_convertible<
	 typename T::independent_variable_type,
	 scalar_trait_t<typename T::state_type> >::value
  >
  > : std::true_type{};

#endif


#ifdef PRESSIO_ENABLE_CXX20
template<class T>
concept FDJ_Common =
  requires { typename T::independent_variable_type; } &&
  requires { typename T::state_type; } &&
  requires { typename T::discrete_residual_type; } &&
  requires { typename T::discrete_jacobian_type; } &&
  std::copy_constructible<typename T::state_type> &&
  std::copy_constructible<typename T::discrete_residual_type> &&
  std::copy_constructible<typename T::discrete_jacobian_type> &&
  requires (const T& A) {
    { A.createState()            } -> std::same_as<typename T::state_type>;
    { A.createDiscreteResidual() } -> std::same_as<typename T::discrete_residual_type>;
    { A.createDiscreteJacobian() } -> std::same_as<typename T::discrete_jacobian_type>;
  };

// ---- Short aliases to keep requires-expressions readable ---------------------
template<class T> using step_t = typename ::pressio::ode::StepCount::value_type; // adjust if you have a different step type
template<class T> using time_t = typename T::independent_variable_type;
template<class T> using state_t = typename T::state_type;
template<class T> using res_t   = typename T::discrete_residual_type;
template<class T> using jac_t   = typename T::discrete_jacobian_type;
template<class T> using opt_jac_ptr_t = std::optional<jac_t<T>*>;

// ---- Method presence checks for NumStates = 1,2,3,4 --------------------------
template<class T>
concept HasDRJ_1 =
  requires (const T& a) {
    { a.discreteResidualAndJacobian(
        std::declval<step_t<T> const&>(),
        std::declval<time_t<T> const&>(),
        std::declval<time_t<T> const&>(),
        std::declval<res_t<T>&>(),
        std::declval<opt_jac_ptr_t<T>>(),
        std::declval<state_t<T> const&>()
      )
    } -> std::same_as<void>;
  };

template<class T>
concept HasDRJ_2 =
  requires (const T& a) {
    { a.discreteResidualAndJacobian(
        std::declval<step_t<T> const&>(),
        std::declval<time_t<T> const&>(),
        std::declval<time_t<T> const&>(),
        std::declval<res_t<T>&>(),
        std::declval<opt_jac_ptr_t<T>>(),
        std::declval<state_t<T> const&>(),
        std::declval<state_t<T> const&>()
      )
    } -> std::same_as<void>;
  };

template<class T>
concept HasDRJ_3 =
  requires (const T& a) {
    { a.discreteResidualAndJacobian(
        std::declval<step_t<T> const&>(),
        std::declval<time_t<T> const&>(),
        std::declval<time_t<T> const&>(),
        std::declval<res_t<T>&>(),
        std::declval<opt_jac_ptr_t<T>>(),
        std::declval<state_t<T> const&>(),
        std::declval<state_t<T> const&>(),
        std::declval<state_t<T> const&>()
      )
    } -> std::same_as<void>;
  };

template<class T>
concept HasDRJ_4 =
  requires (const T& a) {
    { a.discreteResidualAndJacobian(
        std::declval<step_t<T> const&>(),
        std::declval<time_t<T> const&>(),
        std::declval<time_t<T> const&>(),
        std::declval<res_t<T>&>(),
        std::declval<opt_jac_ptr_t<T>>(),
        std::declval<state_t<T> const&>(),
        std::declval<state_t<T> const&>(),
        std::declval<state_t<T> const&>(),
        std::declval<state_t<T> const&>()
      )
    } -> std::same_as<void>;
  };

// select the right arity by NumStates
template<class T, int NumStates>
concept HasDRJ_N =
  (NumStates == 1 && HasDRJ_1<T>) ||
  (NumStates == 2 && HasDRJ_2<T>) ||
  (NumStates == 3 && HasDRJ_3<T>) ||
  (NumStates == 4 && HasDRJ_4<T>);

// Actual concepts
template <class T, int NumStates>
concept FullyDiscreteSystemWithJacobian =
  FDJ_Common<T> && HasDRJ_N<T, NumStates>;

template <class T, int NumStates>
concept RealValuedFullyDiscreteSystemWithJacobian =
     FullyDiscreteSystemWithJacobian<T, NumStates>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::discrete_residual_type> >
  && std::floating_point< scalar_trait_t<typename T::discrete_jacobian_type> >
  && std::convertible_to<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >;

#else

template<class T, int NumStates,class enable = void>
struct FullyDiscreteSystemWithJacobian : std::false_type{};

template<class T, int NumStates>
struct FullyDiscreteSystemWithJacobian<
  T, NumStates,
  std::enable_if_t<
       ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_discrete_residual_typedef<T>::value
    && ::pressio::has_discrete_jacobian_typedef<T>::value
    && std::is_copy_constructible<typename T::state_type>::value
    && std::is_copy_constructible<typename T::discrete_residual_type>::value
    && std::is_copy_constructible<typename T::discrete_jacobian_type>::value
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type >::value
    && ::pressio::ode::has_const_create_discrete_residual_method_return_result<
      T, typename T::discrete_residual_type>::value
    && ::pressio::ode::has_const_create_discrete_jacobian_method_return_result<
      T, typename T::discrete_jacobian_type>::value
    //
    && ::pressio::ode::has_const_discrete_residual_jacobian_method<
      T, NumStates,
      typename ::pressio::ode::StepCount::value_type,
      typename T::independent_variable_type,
      typename T::state_type,
      typename T::discrete_residual_type,
      typename T::discrete_jacobian_type
      >::value
    >
  > : std::true_type{};

template <class T, int NumStates, class = void>
struct RealValuedFullyDiscreteSystemWithJacobian : std::false_type{};

template <class T, int NumStates>
struct RealValuedFullyDiscreteSystemWithJacobian<
  T, NumStates,
  std::enable_if_t<
    FullyDiscreteSystemWithJacobian<T, NumStates>::value
    && std::is_floating_point< scalar_trait_t<typename T::state_type> >::value
    && std::is_floating_point< scalar_trait_t<typename T::discrete_residual_type> >::value
    && std::is_floating_point< scalar_trait_t<typename T::discrete_jacobian_type> >::value
    && std::is_convertible<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >::value
  >
  > : std::true_type{};

#endif






//
// ImplicitResidualJacobianPolicy
//
#ifdef PRESSIO_ENABLE_CXX20

template<class T>
concept ImplicitResidualJacobianPolicy =
  // --- required nested types
  requires {
    typename T::independent_variable_type;
    typename T::state_type;
    typename T::residual_type;
    typename T::jacobian_type;
  }
  // --- known pressio-ops data types & scalar compatibility
  && ::pressio::ops::is_known_data_type<typename T::state_type>::value
  && ::pressio::ops::is_known_data_type<typename T::residual_type>::value
  && ::pressio::ops::is_known_data_type<typename T::jacobian_type>::value
  && ::pressio::all_have_traits_and_same_scalar<
       typename T::state_type,
       typename T::residual_type,
       typename T::jacobian_type
     >::value
  && std::convertible_to<
       typename T::independent_variable_type,
       scalar_trait_t<typename T::state_type>
     >
  // --- factories on a const object
  && requires (T const& a) {
       { a.createState()    } -> std::same_as<typename T::state_type>;
       { a.createResidual() } -> std::same_as<typename T::residual_type>;
       { a.createJacobian() } -> std::same_as<typename T::jacobian_type>;
     }
  // --- call operator with exact signature; must return void and be const
  && requires (T const& a) {
       { a(
           std::declval<StepScheme const&>(),
           std::declval<typename T::state_type const&>(),
           std::declval<ImplicitStencilStatesDynamicContainer<typename T::state_type> const&>(),
           std::declval<ImplicitStencilRightHandSideDynamicContainer<typename T::residual_type>&>(),
           std::declval<::pressio::ode::StepEndAt<typename T::independent_variable_type>>(),
           std::declval<::pressio::ode::StepCount>(),
           std::declval<::pressio::ode::StepSize<typename T::independent_variable_type>>(),
           std::declval<typename T::residual_type &>(),
           std::declval<std::optional<typename T::jacobian_type*> >()
         ) } -> std::same_as<void>;
     };

#else

template<class T, class = void>
struct ImplicitResidualJacobianPolicy : std::false_type{};

template<class T>
struct ImplicitResidualJacobianPolicy<
  T,
  std::enable_if_t<
    ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && ::pressio::has_residual_typedef<T>::value
    && ::pressio::has_jacobian_typedef<T>::value
    //
    && ::pressio::ops::is_known_data_type<typename T::state_type>::value
    && ::pressio::ops::is_known_data_type<typename T::residual_type>::value
    && ::pressio::ops::is_known_data_type<typename T::jacobian_type>::value
    && all_have_traits_and_same_scalar<
      typename T::state_type,
      typename T::residual_type,
      typename T::jacobian_type>::value
    && std::is_convertible<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type>>::value
    //
    // create methods
    //
    && ::pressio::ode::has_const_create_state_method_return_result<
      T, typename T::state_type>::value
    && std::is_same<
      typename T::residual_type,
      decltype(std::declval<T const>().createResidual())
      >::value
    && std::is_same<
      typename T::jacobian_type,
      decltype(std::declval<T const>().createJacobian())
      >::value

    && std::is_void<
      decltype
      (
       std::declval<T const>()
       (
	std::declval<StepScheme const &>(),
	std::declval<typename T::state_type const &>(),
	std::declval<ImplicitStencilStatesDynamicContainer<typename T::state_type> const & >(),
	std::declval<ImplicitStencilRightHandSideDynamicContainer<typename T::residual_type> & >(),
	std::declval< ::pressio::ode::StepEndAt<typename T::independent_variable_type> >(),
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepSize<typename T::independent_variable_type> >(),
	std::declval<typename T::residual_type &>(),
	std::declval< std::optional<typename T::jacobian_type*> >()
	)
       )
      >::value
    >
  > : std::true_type{};

#endif // end ifdef CXX20



#ifdef PRESSIO_ENABLE_CXX20

template<class T>
concept StepperWithoutSolver =
  // required nested types
  requires {
    typename T::independent_variable_type;
    typename T::state_type;
  }
  // canonical call must exist and return void
  && requires(
       T& stepper,
       typename T::state_type& y,
       const ::pressio::ode::StepStartAt<typename T::independent_variable_type>& t0,
       const ::pressio::ode::StepCount& k,
       const ::pressio::ode::StepSize<typename T::independent_variable_type>& dt
     ) {
       { stepper(y, t0, k, dt) } -> std::same_as<void>;
     }
  // forbid accepting const state&
  && (!requires(
        T& stepper,
        const typename T::state_type& y,
        const ::pressio::ode::StepStartAt<typename T::independent_variable_type>& t0,
        const ::pressio::ode::StepCount& k,
        const ::pressio::ode::StepSize<typename T::independent_variable_type>& dt
      ) {
        stepper(y, t0, k, dt);
      })
  // forbid accepting state&&
  && (!requires(
        T& stepper,
        typename T::state_type&& y,
        const ::pressio::ode::StepStartAt<typename T::independent_variable_type>& t0,
        const ::pressio::ode::StepCount& k,
        const ::pressio::ode::StepSize<typename T::independent_variable_type>& dt
      ) {
        stepper(std::move(y), t0, k, dt);
      });

#else

/*
  This comment documents the pieces defined just below:

    1) `stepper_call_t<T>` — an alias for the *deduced return type* of calling
       `T::operator()` with the expected arguments.
    2) `StepperWithoutSolver<T>` — a boolean type trait that evaluates to
       `true` iff that call is well-formed and returns **exactly `void`**.

  ---------------------------------------------------------------------------
  How does it work?
  ---------------------------------------------------------------------------
  1) `stepper_call_t<T>` forms the type of the expression

       std::declval<T&>()(
         std::declval<typename T::state_type&>(),
         std::declval<StepStartAt<typename T::independent_variable_type>>(),
         std::declval<StepCount>(),
         std::declval<StepSize<typename T::independent_variable_type>>()
       )

     using `decltype(...)`. `std::declval<X>()` fabricates a value of type X in an
     unevaluated context, letting us ask “what would this call return?”. If the call
     is ill-formed (wrong/missing types, or requires non-const lvalue refs for the
     last three arguments), `stepper_call_t<T>` itself fails to form and the trait
     is SFINAE’d out.

  2) The partial specialization of `StepperWithoutSolver<T, ...>` is gated by
     `std::void_t<...>`:

       std::void_t<
         typename T::independent_variable_type,
         typename T::state_type,
         stepper_call_t<T>
       >

     This specialization only participates if the nested typedefs exist *and*
     the call expression type could be formed.

  3) Inside the participating specialization we apply additional semantic checks
     with a `std::conjunction`:

     - `std::is_same<void, stepper_call_t<T>>` requires the deduced return type
       to be **exactly `void`**.

     - A positive invocability check:
         std::is_invocable_r<
           void, T&,
           state_type&,
           const StepStartAt<IV>&, const StepCount&, const StepSize<IV>&>
       This enforces that:
         * `state` is accepted as a **non-const lvalue reference** (so it can be
           modified).
         * The other three parameters are accepted as **non-modifiable** inputs.
           Requiring `const&` here also permits implementations that take them
           **by value** (since prvalues bind to `const&`).

     - Two negative invocability checks:
         * Not invocable with `const state_type&` (rejects implementations that
           cannot modify `state`).
         * Not invocable with `state_type&&` (rejects by-value / forwarding-ref
           signatures for `state`, which could allow moving/copying instead of
           in-place modification).
*/

template<class T>
using stepper_call_t = decltype(
  std::declval<T&>()(
    std::declval<typename T::state_type &>(),
    std::declval<StepStartAt<typename T::independent_variable_type>>(),
    std::declval<StepCount>(),
    std::declval<StepSize<typename T::independent_variable_type>>()
  )
);

template <class T, class = void>
struct StepperWithoutSolver : std::false_type{};

template <class T>
struct StepperWithoutSolver<
  T,
  std::void_t<
    typename T::independent_variable_type,
    typename T::state_type,
    stepper_call_t<T>
    >
> : std::conjunction<
    // must return void for the canonical call
    std::is_same<void, stepper_call_t<T>>,

    // state MUST be a non-const lvalue reference (callable with state&)
    std::is_invocable_r<
      void, T&,
      typename T::state_type&,
      // other params must be const-ref OR by-value: requiring const& works for both
      const StepStartAt<typename T::independent_variable_type>&,
      const StepCount&,
      const StepSize<typename T::independent_variable_type>&
    >,

    // forbid accepting a const state& (would not allow modification)
    std::negation<std::is_invocable<
      T&,
      const typename T::state_type&,
      const StepStartAt<typename T::independent_variable_type>&,
      const StepCount&,
      const StepSize<typename T::independent_variable_type>&
    >>,

    // forbid accepting state&& (would allow by-value / move-only patterns)
    std::negation<std::is_invocable<
      T&,
      typename T::state_type&&,
      const StepStartAt<typename T::independent_variable_type>&,
      const StepCount&,
      const StepSize<typename T::independent_variable_type>&
    >>
  > {};

#endif


/* -------------------------------------------------------------
  StateObserver
*/
#ifdef PRESSIO_ENABLE_CXX20

template<class T, class IndVarType, class StateType>
concept StateObserver =
  requires (const T& obs, StepCount k, IndVarType t, const StateType& x)
  {
    { obs(k, t, x) } -> std::same_as<void>;
  };

#else

/*
  Trait to detect a **const** state observer callable.

  Succeeds iff a type `T` has a **const** call operator with the shape:

    void operator()( ::pressio::ode::StepCount,
                     IndVarType,
                     StateType const& ) const;

  More precisely:
   - We form the type of the expression
       std::declval<T const&>()(StepCount{}, IndVarType{}, StateType const&{})
     via `decltype(...)`. If that call is ill-formed, the specialization
     is discarded (SFINAE) and the trait is `false`.
   - If the call is well-formed, we additionally require the deduced
     return type to be **exactly `void`** (not just convertible to void).

  Usage:
    static_assert(StateObserver<MyObs, double, MyState>::value, "Must be a const observer");
*/

template<class T, class IndVarType, class StateType>
using state_observer_call_c_t = decltype(
  std::declval<T const&>()(
    std::declval<StepCount>(),
    std::declval<IndVarType>(),
    std::declval<StateType const&>()
  )
);

template<class T, class IndVarType, class StateType, class = void>
struct StateObserver : std::false_type {};

template<class T, class IndVarType, class StateType>
struct StateObserver<
  T, IndVarType, StateType,
  std::void_t< state_observer_call_c_t<T, IndVarType, StateType> >
>
: std::is_same<void, state_observer_call_c_t<T, IndVarType, StateType>> {};

#endif


/* -------------------------------------------------------------
   RhsObserver
*/
#ifdef PRESSIO_ENABLE_CXX20

template<class T, class IndVarType, class RhsType>
concept RhsObserver =
  requires (const T& obs,
            StepCount k,
            IntermediateStepCount subk,
            IndVarType t,
            const RhsType& r)
  {
    { obs(k, subk, t, r) } -> std::same_as<void>;
  };

#else

/*
   Detect a **const** RHS observer with exact-void return:
   void operator()(StepCount, IntermediateStepCount, IndVarType, RhsType const&) const;
*/

template<class T, class IndVarType, class RhsType>
using rhs_observer_call_c_t = decltype(
  std::declval<T const&>()(
    std::declval<StepCount>(),
    std::declval<IntermediateStepCount>(),
    std::declval<IndVarType>(),
    std::declval<RhsType const&>()
  )
);

// Primary: default to false
template<class T, class IndVarType, class RhsType, class = void>
struct RhsObserver : std::false_type {};

// Specialization: only if the const-call is well-formed; require return type == void
template<class T, class IndVarType, class RhsType>
struct RhsObserver<
  T, IndVarType, RhsType,
  std::void_t< rhs_observer_call_c_t<T, IndVarType, RhsType> >
>
: std::is_same<void, rhs_observer_call_c_t<T, IndVarType, RhsType>> {};

#endif


/* -------------------------------------------------------------
  StepSizePolicy
*/
#ifdef PRESSIO_ENABLE_CXX20

template<class T, class IndVarType>
concept StepSizePolicy =
  // must be callable with non-const lvalue ref and return void
  requires (const T& policy,
            StepCount k,
            StepStartAt<IndVarType> start,
            StepSize<IndVarType>& dt) {
    { policy(k, start, dt) } -> std::same_as<void>;
  }
  // must NOT be callable with const lvalue ref
  && (!requires (const T& policy,
                 StepCount k,
                 StepStartAt<IndVarType> start) {
        policy(k, start,
               std::declval<const StepSize<IndVarType>&>());
     })
  // must NOT be callable with rvalue
  && (!requires (const T& policy,
                 StepCount k,
                 StepStartAt<IndVarType> start) {
        policy(k, start,
               std::declval<StepSize<IndVarType>&&>());
     });

#else

/*
  Detect at compile time whether a type `T` can act as a *step size policy*,
  i.e., it is a callable object with a **const** call operator of the form:

      void operator()( ::pressio::ode::StepCount,
                       ::pressio::ode::StepStartAt<IndVarType>,
                       ::pressio::ode::StepSize<IndVarType>& ) const;

  Semantics
  ---------
  - The third argument is a non-const reference to StepSize<IndVarType>, so the
    policy can write/adjust the step size.
  - The trait requires the return type to be **exactly `void`**.
*/

template<class T, class IndVarType>
using step_size_policy_call_c_t = decltype(
  std::declval<T const&>()(
    std::declval<StepCount>(),
    std::declval<StepStartAt<IndVarType>>(),
    std::declval<StepSize<IndVarType>&>()
  )
);

template <class T, class IndVarType, class Enable = void>
struct StepSizePolicy : std::false_type {};

template <class T, class IndVarType>
struct StepSizePolicy<
  T, IndVarType,
  std::void_t< step_size_policy_call_c_t<T, IndVarType> >
>
: std::conjunction<
    // must be callable and return void with non-const lvalue ref
    std::is_same<void, step_size_policy_call_c_t<T, IndVarType>>,
    std::is_invocable_r<
      void, const T&, StepCount, StepStartAt<IndVarType>, StepSize<IndVarType>&>,
    // must NOT be callable with const lvalue ref
    std::negation<std::is_invocable<
      const T&, StepCount, StepStartAt<IndVarType>, const StepSize<IndVarType>&>>,
    // must NOT be callable with rvalue
    std::negation<std::is_invocable<
      const T&, StepCount, StepStartAt<IndVarType>, StepSize<IndVarType>&&>>
  > {};

#endif



/* -------------------------------------------------------------
  StepSizePolicyWithReductionScheme
*/
#ifdef PRESSIO_ENABLE_CXX20

template<class T, class IndVarType>
concept StepSizePolicyWithReductionScheme =
  // must be callable and return void with non-const lvalue refs
  requires (const T& policy,
            StepCount k,
            StepStartAt<IndVarType> start,
            StepSize<IndVarType>& dt,
            StepSizeMinAllowedValue<IndVarType>& dt_min,
            StepSizeScalingFactor<IndVarType>& scale) {
    { policy(k, start, dt, dt_min, scale) } -> std::same_as<void>;
  }
  // must NOT accept any of the three by const&
  && (!requires (const T& policy,
                 StepCount k,
                 StepStartAt<IndVarType> start,
                 const StepSize<IndVarType>& dtc,
                 StepSizeMinAllowedValue<IndVarType>& dt_min,
                 StepSizeScalingFactor<IndVarType>& scale) {
        policy(k, start, dtc, dt_min, scale);
     })
  && (!requires (const T& policy,
                 StepCount k,
                 StepStartAt<IndVarType> start,
                 StepSize<IndVarType>& dt,
                 const StepSizeMinAllowedValue<IndVarType>& dt_min_c,
                 StepSizeScalingFactor<IndVarType>& scale) {
        policy(k, start, dt, dt_min_c, scale);
     })
  && (!requires (const T& policy,
                 StepCount k,
                 StepStartAt<IndVarType> start,
                 StepSize<IndVarType>& dt,
                 StepSizeMinAllowedValue<IndVarType>& dt_min,
                 const StepSizeScalingFactor<IndVarType>& scale_c) {
        policy(k, start, dt, dt_min, scale_c);
     })
  // must NOT accept any of the three by rvalue (disallow by-value)
  && (!requires (const T& policy,
                 StepCount k,
                 StepStartAt<IndVarType> start,
                 StepSize<IndVarType>&& dtr,
                 StepSizeMinAllowedValue<IndVarType>& dt_min,
                 StepSizeScalingFactor<IndVarType>& scale) {
        policy(k, start, std::move(dtr), dt_min, scale);
     })
  && (!requires (const T& policy,
                 StepCount k,
                 StepStartAt<IndVarType> start,
                 StepSize<IndVarType>& dt,
                 StepSizeMinAllowedValue<IndVarType>&& dt_min_r,
                 StepSizeScalingFactor<IndVarType>& scale) {
        policy(k, start, dt, std::move(dt_min_r), scale);
     })
  && (!requires (const T& policy,
                 StepCount k,
                 StepStartAt<IndVarType> start,
                 StepSize<IndVarType>& dt,
                 StepSizeMinAllowedValue<IndVarType>& dt_min,
                 StepSizeScalingFactor<IndVarType>&& scale_r) {
        policy(k, start, dt, dt_min, std::move(scale_r));
     });

#else

/*
  Detect at compile time whether a type `T` can act as a *step size policy with
  reduction scheme*, i.e., it has a **const** call operator with the shape:

      void operator()( ::pressio::ode::StepCount,
                       ::pressio::ode::StepStartAt<IndVarType>,
                       ::pressio::ode::StepSize<IndVarType>&,
                       ::pressio::ode::StepSizeMinAllowedValue<IndVarType>&,
                       ::pressio::ode::StepSizeScalingFactor<IndVarType>& ) const;

  Semantics
  ---------
  - The last three parameters are non-const references the policy can modify:
      * StepSize<IndVarType>&                — current step size
      * StepSizeMinAllowedValue<IndVarType>& — minimum allowed step size
      * StepSizeScalingFactor<IndVarType>&   — multiplicative scaling factor
  - The trait requires the return type to be **exactly `void`**.
*/

template<class T, class IndVarType>
using step_size_policy_with_reduction_call_c_t = decltype(
  std::declval<T const&>()(
    std::declval<StepCount>(),
    std::declval<StepStartAt<IndVarType>>(),
    std::declval<StepSize<IndVarType>&>(),
    std::declval<StepSizeMinAllowedValue<IndVarType>&>(),
    std::declval<StepSizeScalingFactor<IndVarType>&>()
  )
);

template <class T, class IndVarType, class Enable = void>
struct StepSizePolicyWithReductionScheme : std::false_type {};

template <class T, class IndVarType>
struct StepSizePolicyWithReductionScheme<
  T, IndVarType,
  std::void_t< step_size_policy_with_reduction_call_c_t<T, IndVarType> >
>
: std::conjunction<
    // must return void and be callable with non-const lvalue refs
    std::is_same<void, step_size_policy_with_reduction_call_c_t<T, IndVarType>>,
    std::is_invocable_r<
      void, const T&,
      StepCount, StepStartAt<IndVarType>,
      StepSize<IndVarType>&,
      StepSizeMinAllowedValue<IndVarType>&,
      StepSizeScalingFactor<IndVarType>&
    >,
    // must NOT be callable with any of the three by const&
    std::negation<std::is_invocable<
      const T&, StepCount, StepStartAt<IndVarType>,
      const StepSize<IndVarType>&,
      StepSizeMinAllowedValue<IndVarType>&,
      StepSizeScalingFactor<IndVarType>&
    >>,
    std::negation<std::is_invocable<
      const T&, StepCount, StepStartAt<IndVarType>,
      StepSize<IndVarType>&,
      const StepSizeMinAllowedValue<IndVarType>&,
      StepSizeScalingFactor<IndVarType>&
    >>,
    std::negation<std::is_invocable<
      const T&, StepCount, StepStartAt<IndVarType>,
      StepSize<IndVarType>&,
      StepSizeMinAllowedValue<IndVarType>&,
      const StepSizeScalingFactor<IndVarType>&
    >>,
    // must NOT be callable with any of the three by rvalue (would allow by-value)
    std::negation<std::is_invocable<
      const T&, StepCount, StepStartAt<IndVarType>,
      StepSize<IndVarType>&&,
      StepSizeMinAllowedValue<IndVarType>&,
      StepSizeScalingFactor<IndVarType>&
    >>,
    std::negation<std::is_invocable<
      const T&, StepCount, StepStartAt<IndVarType>,
      StepSize<IndVarType>&,
      StepSizeMinAllowedValue<IndVarType>&&,
      StepSizeScalingFactor<IndVarType>&
    >>,
    std::negation<std::is_invocable<
      const T&, StepCount, StepStartAt<IndVarType>,
      StepSize<IndVarType>&,
      StepSizeMinAllowedValue<IndVarType>&,
      StepSizeScalingFactor<IndVarType>&&
    >>
  > {};

#endif

}}
#endif  // PRESSIOROM_ODE_ODE_CONCEPTS_HPP_
