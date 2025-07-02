/*
//@HEADER
// ************************************************************************
//
// solvers_nonlinear_enums.hpp
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

#ifndef PRESSIOROM_SOLVERS_NONLINEAR_SOLVERS_NONLINEAR_ENUMS_AND_TAGS_HPP_
#define PRESSIOROM_SOLVERS_NONLINEAR_SOLVERS_NONLINEAR_ENUMS_AND_TAGS_HPP_

namespace pressio{ namespace nonlinearsolvers{

enum class Stop{
  WhenAbsolutel2NormOfCorrectionBelowTolerance,
  WhenRelativel2NormOfCorrectionBelowTolerance,
  WhenAbsolutel2NormOfResidualBelowTolerance,
  WhenRelativel2NormOfResidualBelowTolerance,
  WhenAbsolutel2NormOfGradientBelowTolerance,
  WhenRelativel2NormOfGradientBelowTolerance,
  WhenAbsoluteObjectiveBelowTolerance,
  WhenRelativeObjectiveBelowTolerance,
  AfterMaxIters
};

enum class Update{
  Standard,
  Armijo,
  BacktrackStrictlyDecreasingObjective,
  LMSchedule1,
  LMSchedule2,
  Custom
};

enum class Diagnostic{
  correctionAbsolutel2Norm,
  correctionRelativel2Norm,
  residualAbsolutel2Norm,
  residualRelativel2Norm,
  gradientAbsolutel2Norm,
  gradientRelativel2Norm,
  objectiveAbsolute,
  objectiveRelative,
  invalid
};

template<class T = void>
std::string diagnostic_to_string(Diagnostic d){
  switch(d)
    {
    case Diagnostic::correctionRelativel2Norm: return "correctionRelativel2Norm";
    case Diagnostic::correctionAbsolutel2Norm: return "correctionAbsolutel2Norm";
    case Diagnostic::residualRelativel2Norm:   return "residualRelativel2Norm";
    case Diagnostic::residualAbsolutel2Norm:   return "residualAbsolutel2Norm";
    case Diagnostic::gradientRelativel2Norm:   return "gradientRelativel2Norm";
    case Diagnostic::gradientAbsolutel2Norm:   return "gradientAbsolutel2Norm";
    case Diagnostic::objectiveAbsolute:	       return "objectiveAbsolute";
    case Diagnostic::objectiveRelative:	       return "objectiveRelative";
    default: return "invalid";
    };
};

template<class T = void>
bool is_absolute_diagnostic(Diagnostic d){
  switch(d)
    {
    case Diagnostic::correctionRelativel2Norm: return false;
    case Diagnostic::correctionAbsolutel2Norm: return true;
    case Diagnostic::residualRelativel2Norm:   return false;
    case Diagnostic::residualAbsolutel2Norm:   return true;
    case Diagnostic::gradientRelativel2Norm:   return false;
    case Diagnostic::gradientAbsolutel2Norm:   return true;
    case Diagnostic::objectiveAbsolute:   return true;
    case Diagnostic::objectiveRelative:   return false;
    default: return true;
    };
};

template<class T = void>
std::string diagnostic_to_log_symbol(Diagnostic d){
  switch(d)
    {
    case Diagnostic::correctionRelativel2Norm: return "||delta||(r)";
    case Diagnostic::correctionAbsolutel2Norm: return "||delta||(a)";
    case Diagnostic::residualRelativel2Norm:   return "||R||(r)";
    case Diagnostic::residualAbsolutel2Norm:   return "||R||(a)";
    case Diagnostic::gradientRelativel2Norm:   return "||g||(r)";
    case Diagnostic::gradientAbsolutel2Norm:   return "||g||(a)";
    case Diagnostic::objectiveAbsolute:   return "obj(a)";
    case Diagnostic::objectiveRelative:   return "obj(r)";
    default:  return "invalid";
    };
};

template<class T = void>
std::string diagnostic_to_log_format(Diagnostic d){
  switch(d)
    {
    case Diagnostic::correctionRelativel2Norm:
    case Diagnostic::correctionAbsolutel2Norm:
    case Diagnostic::residualRelativel2Norm:
    case Diagnostic::residualAbsolutel2Norm:
    case Diagnostic::gradientRelativel2Norm:
    case Diagnostic::gradientAbsolutel2Norm:
    case Diagnostic::objectiveAbsolute:
    case Diagnostic::objectiveRelative:
      return "{:.6e}";

    default:
      return "invalid";
    };
};

template<class T = void>
Diagnostic stop_criterion_to_public_diagnostic(const Stop & sc)
{
  switch (sc)
    {
    case Stop::WhenAbsolutel2NormOfCorrectionBelowTolerance: return Diagnostic::correctionAbsolutel2Norm;
    case Stop::WhenRelativel2NormOfCorrectionBelowTolerance: return Diagnostic::correctionRelativel2Norm;
    case Stop::WhenAbsolutel2NormOfResidualBelowTolerance: return Diagnostic::residualAbsolutel2Norm;
    case Stop::WhenRelativel2NormOfResidualBelowTolerance: return Diagnostic::residualRelativel2Norm;
    case Stop::WhenAbsolutel2NormOfGradientBelowTolerance: return Diagnostic::gradientAbsolutel2Norm;
    case Stop::WhenRelativel2NormOfGradientBelowTolerance: return Diagnostic::gradientRelativel2Norm;
    case Stop::WhenAbsoluteObjectiveBelowTolerance: return Diagnostic::objectiveAbsolute;
    case Stop::WhenRelativeObjectiveBelowTolerance: return Diagnostic::objectiveRelative;
    default: return Diagnostic::invalid;
    };
}

/*
  List of tags used to label data in
  the solver data registry
*/

// the state is the "solution", what is updated
// and what we solve for
struct StateTag{};

// the solver for the newton step
struct InnerSolverTag{};

// any newton solve is an iterative process where the
// next state x_k+1 is computed as: x_k+1 = x_k + alpha*delta
// where delta is the correction. Here, the correction
// is therefore defined as: delta = x_k+1 - x_k
struct CorrectionTag{};

// self-explanatory
struct ResidualTag{};
struct JacobianTag{};
struct GradientTag{};
struct HessianTag{};

// a "scratch" state that can be used as trial state
// needed when doing a line search
struct LineSearchTrialStateTag{};

// the damping factor for LM
struct LevenbergMarquardtDampingTag{};
// for LM, we need to store the unscaled hessian
struct LevenbergMarquardtUndampedHessianTag{};

// when doing, e.g., a weighted least-squares, we need
// access to the weighting operator
struct WeightingOperatorTag{};

// self-explanatory
struct WeightedResidualTag{};
struct WeightedJacobianTag{};

}} //end namespace pressio::nonlinearsolvers
#endif  // PRESSIOROM_SOLVERS_NONLINEAR_SOLVERS_NONLINEAR_ENUMS_AND_TAGS_HPP_
