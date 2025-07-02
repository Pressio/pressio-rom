/*
//@HEADER
// ************************************************************************
//
// ode_enum_and_tags.hpp
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

#ifndef PRESSIOROM_ODE_ODE_ENUM_AND_TAGS_HPP_
#define PRESSIOROM_ODE_ODE_ENUM_AND_TAGS_HPP_

namespace pressio{ namespace ode{

/**
 * Enum class listing all supported ODE step schemes.
 *
 * Explicit schemes:
 *   - ForwardEuler: First-order explicit Euler scheme
 *   - RungeKutta4: Classical fourth-order Runge-Kutta scheme
 *   - AdamsBashforth2: Second-order Adams-Bashforth multi-step method
 *   - SSPRungeKutta3: Strong Stability Preserving third-order Runge-Kutta
 *
 * Implicit schemes:
 *   - BDF1: First-order Backward Differentiation Formula
 *   - BDF2: Second-order Backward Differentiation Formula
 *   - CrankNicolson: Second-order Crank-Nicolson scheme
 *   - ImplicitArbitrary: Placeholder for arbitrary user-defined implicit schemes
 */

enum class StepScheme{
  // explicit
  ForwardEuler,
  RungeKutta4,
  AdamsBashforth2,
  SSPRungeKutta3,
  // implicit
  BDF1,
  BDF2,
  CrankNicolson,
  ImplicitArbitrary
};

/**
 * Utility: Checks if a StepScheme is explicit.
 * Returns true if the scheme is one of the supported explicit schemes.
 */
template<class T = bool>
T is_explicit_scheme(StepScheme name)
{
  if (name == StepScheme::ForwardEuler){ return true; }
  else if (name == StepScheme::RungeKutta4){ return true; }
  else if (name == StepScheme::AdamsBashforth2){ return true; }
  else if (name == StepScheme::SSPRungeKutta3){ return true; }
  else{ return false; }
}

/**
 * Utility: Checks if a StepScheme is implicit.
 * Returns true if the scheme is not explicit.
 */
template<class T = bool>
T is_implicit_scheme(StepScheme name){
  return !is_explicit_scheme(name);
}

/**
 * Tags for time integration schemes.
 *
 * These tags are used both internally for overloading
 * as well as used by other modules, like ROMs.
 */
struct ForwardEuler{};
struct RungeKutta4{};
struct AdamsBashforth2{};
struct SSPRungeKutta3{};
struct BDF1{};
struct BDF2{};
struct CrankNicolson{};
struct ImplicitArbitrary{};

/**
 * Time-level tags for multi-step or multi-stage schemes.
 *
 *   - nPlusOne: Represents the (n+1)-th time level.
 *   - n:        Represents the n-th time level.
 *   - nMinusOne:   Represents the (n-1)-th time level.
 *   - nMinusTwo:   Represents the (n-2)-th time level.
 *   - nMinusThree: Represents the (n-3)-th time level.
 *   - nMinusFour:  Represents the (n-4)-th time level.
 */
class nPlusOne{};
class n{};
class nMinusOne{};
class nMinusTwo{};
class nMinusThree{};

}}//end namespace pressio::ode

#endif  // PRESSIOROM_ODE_ODE_ENUM_AND_TAGS_HPP_
