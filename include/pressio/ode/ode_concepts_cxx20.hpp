/*
//@HEADER
// ************************************************************************
//
// ode_concepts_cxx20.hpp
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

#ifndef PRESSIOROM_ODE_CONCEPTS_ODE_CONCEPTS_CXX20_HPP_
#define PRESSIOROM_ODE_CONCEPTS_ODE_CONCEPTS_CXX20_HPP_

#include <concepts>
#include <optional>

namespace pressio{ namespace ode{

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

//
// refine for real-valued case
//
template <class T>
concept RealValuedOdeSystem =
     OdeSystem<T>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::rhs_type> >
  && std::convertible_to<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >;

template <class T>
concept RealValuedOdeSystemFusingMassMatrixAndRhs =
  OdeSystemFusingMassMatrixAndRhs<T>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::rhs_type> >
  && std::floating_point< scalar_trait_t<typename T::mass_matrix_type> >
  && std::convertible_to<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >;

template <class T>
concept RealValuedOdeSystemFusingRhsAndJacobian =
     OdeSystemFusingRhsAndJacobian<T>
  && std::floating_point< scalar_trait_t<typename T::state_type> >
  && std::floating_point< scalar_trait_t<typename T::rhs_type> >
  && std::floating_point< scalar_trait_t<typename T::jacobian_type> >
  && std::convertible_to<
      typename T::independent_variable_type,
      scalar_trait_t<typename T::state_type> >;

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


//
// RealValuedFullyDiscreteSystemWithJacobian
//

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

// Final concepts
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

}}
#endif  // PRESSIOROM_ODE_CONCEPTS_ODE_SYSTEM_HPP_
