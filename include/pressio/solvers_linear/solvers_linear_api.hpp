/*
//@HEADER
// ************************************************************************
//
// solvers_linear_solver.hpp
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

#ifndef PRESSIO_SOLVERS_LINEAR_SOLVERS_PUBLIC_API_HPP_
#define PRESSIO_SOLVERS_LINEAR_SOLVERS_PUBLIC_API_HPP_

// the following tags are defined here so that they are visible
// by the impl headers included below.
// These tags that can be used to specify the desired
// solver type when instantiating a solver.

namespace pressio{ namespace linearsolvers{

namespace iterative{
struct CG {};
struct LSCG {};
struct Bicgstab {};
struct GMRES{};
}

namespace direct{
struct HouseholderQR {};
struct ColPivHouseholderQR {};
struct PartialPivLU {};
}

}} // end namespace pressio::linearsolver

#include "./impl/solvers_linear_traits.hpp"
#include "./impl/solvers_linear_solver_selector_impl.hpp"

/*
 * The following `impl::Selector<TagType, MatrixType, Args...>` is a compile-time
 * dispatch mechanism that maps a solver tag and matrix type to an impl class.
 * These implementation classes are **thin wrappers** around existing solvers
 * from supported backends such as Eigen, Kokkos.
 * These implementation classes are organized per-backend:
 *   - For Eigen-backed matrices, selection points to impl classes in:
 *       impl/solvers_linear_eigen_{iterative, direct}.hpp
 *   - For Kokkos-backed matrices, the selection would go to different
 *     backend-specific files, if available.
 *
 * The wrappers:
 *   - Expose a uniform API (e.g., `solve(A, b, x)`, `setTolerance(...)`)
 *   - Internally use traits and tags to determine the correct
 *     backend-specific solver objects (e.g., `Eigen::ConjugateGradient`)
 *   - Are enabled only if the corresponding TPL is enabled
 *
 * -----------------------------------------
 * Dispatch Example
 * -----------------------------------------
 * If you pass:
 *   - `pressio::linearsolvers::iterative::CG` as the solver tag
 *   - `Eigen::SparseMatrix<double>` as the matrix type
 *
 * Then the selected implementation will be: impl::EigenIterativeWrapper
 * which interally would wrap and use Eigen::ConjugateGradient.
 *
 * -----------------------------------------
 * Usage example:
 * -----------------------------------------
 *   using MyMatrixType = ...;
 *   pressio::linearsolvers::Solver<
 *     pressio::linearsolvers::iterative::CG, MyMatrixType
 *   > mySolver;
 *
 * ----------------------------------------------------------
 * NOTE: Implementation Types Are Internal
 * ----------------------------------------------------------
 *
 * Users should **not** depend on or attempt to access the concrete type
 * of the selected solver implementation.
 * Users should interact only with the `Solver<TagType, MatrixType, ...>` alias,
 * which abstracts away the backend and implementation details.
 * Treat `Solver<...>` as an opaque type that meets an interface (see below).
 * Do not rely on the name, namespace, or structure of
 * the underlying implementation class.
 *
 * --------------------------
 * Common API of Implementation Classes
 * --------------------------
 *
 * All implementation classes (i.e., the ones selected via the `Solver<Tag, Matrix, ...>`)
 * expose a **uniform API**, designed to be backend-independent and consistent.
 *
 * Every solver wrapper provides:
 *
 *   `void solve(const MatrixType & A, const VectorType & b, VectorType & x) const`
 *      - Solves the linear system A * x = b and stores the result in `x`.
 *      - Matrix and vector types must be backend-compatible: e.g., you cannot
 *        pass kokkos views if you specified an Eigen-based solver.
 *
 * Additional methods differ between iterative and direct solvers:
 *
 * ---------------------------------
 * API: Iterative Solver Wrappers
 * ---------------------------------
 * These methods are typically available:
 *
 *  iteration_type numIterationsExecuted() const
 *      - get the number of iterations executed
 *
 *  iteration_type maxIterations() const const
 *       - get the maximum number of iterations.
 *
 *  void setMaxIterations(iteration_type maxIters)
 *       - set the maximum number of iterations.
 *
 * ---------------------------------
 * API: Direct Solver Wrappers
 * ---------------------------------
 * These solvers generally expose only:
 *
 *   `void solve(const MatrixType & A, const VectorType & b, VectorType & x) const`
 *
 * They **do not** expose iteration or tolerances, since these are not applicable.
 *
 */
namespace pressio{ namespace linearsolvers{

template<typename TagType, typename MatrixType, typename ... Args>
using Solver = typename impl::Selector<TagType, MatrixType, Args...>::type;

}}
#endif
