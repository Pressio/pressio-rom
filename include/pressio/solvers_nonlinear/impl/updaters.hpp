/*
//@HEADER
// ************************************************************************
//
// updaters.hpp
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

#ifndef PRESSIOROM_SOLVERS_NONLINEAR_IMPL_UPDATERS_HPP_
#define PRESSIOROM_SOLVERS_NONLINEAR_IMPL_UPDATERS_HPP_

namespace pressio{
namespace nlsol{
namespace impl{

class DefaultUpdater{
public:
  template<class RegType, class ... Args>
  void operator()(RegType & reg, Args && ... /*unused*/)
  {
    PRESSIOLOG_DEBUG("nonlinsolver: default update");
    auto & state = reg.template get<StateTag>();
    const auto & c = reg.template get<CorrectionTag>();
    ::pressio::ops::update(state, 1, c, 1);
  }
};

class ArmijoUpdater
{
public:
  template<class RegistryType, class ObjF, class ScalarType>
  void operator()(RegistryType & reg,
		  ObjF objective,
		  ScalarType objectiveValueAtCurrentNewtonStep)
  {
    using scalar_type = std::remove_const_t<ScalarType>;
    PRESSIOLOG_DEBUG("Armijo update");

    // https://people.maths.ox.ac.uk/hauser/hauser_lecture2.pdf

    // Armijo rule says to backtrack alpha until:
    // f(x_trial) - f(x_k) <= alpha_l * beta * dot(g_k, p_k)
    //
    // where:
    // k = the GN step
    // l = indexes the Armijo backtracking stages
    // p_k is the correction at GN k-th step
    // g_k is the gradient at GN k-th step
    // x_trial = x_k + alpha_l * p_k

    const auto & p_k   = reg.template get<CorrectionTag>();
    const auto & g_k   = reg.template get<GradientTag>();
    auto & state = reg.template get<StateTag>();
    auto & x_trial  = reg.template get<LineSearchTrialStateTag>();

    auto alpha = static_cast<scalar_type>(1);
    constexpr auto alpha_lower_bound = static_cast<scalar_type>(0.001);
    const auto beta = static_cast<scalar_type>(0.0001);
    const auto gkDotpk = ::pressio::ops::dot(g_k, p_k);

    PRESSIOLOG_DEBUG("start backtracking");
    constexpr auto zero = static_cast<scalar_type>(0);
    constexpr auto one = static_cast<scalar_type>(1);
    scalar_type ftrial = {};
    while (true)
      {
	if (std::abs(alpha) <= alpha_lower_bound){
	  const auto msg = ": in armijo, alpha = " + std::to_string(alpha)
	    + " <= " + std::to_string(alpha_lower_bound)
	    + ", too small, exiting line search";
	  throw ::pressio::eh::LineSearchStepTooSmall(msg);
	}

	// update : x_trial = x_k + alpha*p_k
	::pressio::ops::update(x_trial, zero, state, one, p_k, alpha);

	// compute rhs_l = alpha_l * beta * dot(g_k, p_k)
	const auto rhs = alpha * beta * gkDotpk;
	const auto lhs = objective(x_trial) - objectiveValueAtCurrentNewtonStep;
	PRESSIOLOG_DEBUG("alpha = {:5f}: (f_trial-f) = {:6e}, rhs = {:6e}", alpha, lhs, rhs);
	if (lhs <= rhs){
	  PRESSIOLOG_DEBUG("condition satisfied: f_trial-f <= rhs, exiting");
	  // solution update: state = state + alpha*p_k
	  ::pressio::ops::update(state, one, p_k, alpha);
	  break;
	}

	// exit when abs(fytrail-fy) < eps, leave eps = 1e-14 for now
	// change later with some machine epsilon
	if (std::abs(lhs) <= 1e-14){
	  const auto msg = ": in armijo, exiting line search";
	  throw ::pressio::eh::LineSearchObjFunctionChangeTooSmall(msg);
	}

	/* convectional way to backtrack:alpha_l+1 = 0.5 * alpha_l */
	alpha *= 0.5;
      }
  }
};

template <typename ScalarType>
class BacktrackStrictlyDecreasingObjectiveUpdater
{
private:
  ScalarType zeta_ = 1;

public:
  BacktrackStrictlyDecreasingObjectiveUpdater(std::optional<std::vector<ScalarType> > parameter)
  {
    if (parameter)
    {
      zeta_ = parameter.value().front();
    }
  }

  template<class RegistryType, class ObjF>
  void operator()(RegistryType & reg,
		  ObjF objective,
		  ScalarType objectiveValueAtCurrentNewtonStep)
  {
    /*
     * Two-phase backtracking line search:
     *   1) Strictly decreasing objective (zeta = 1.0)
     *   2) Relaxed decreasing objective (zeta defined by user)
     *
     * If Phase 1 (strict) fails, we revert to the previous
     * trial state and alpha, and attempt Phase 2 (relaxed).
     */
    using scalar_type = std::remove_const_t<ScalarType>;
    PRESSIOLOG_DEBUG("BacktrackStrictlyDecreasingObjective update");

    auto & trialState  = reg.template get<LineSearchTrialStateTag>();
    auto & state = reg.template get<StateTag>();
    const auto & p_k   = reg.template get<CorrectionTag>();

    auto alpha = static_cast<scalar_type>(1);
    constexpr auto alphaLowerBound = static_cast<scalar_type>(0.001);
    constexpr auto one  = static_cast<scalar_type>(1);
    constexpr auto zero = static_cast<scalar_type>(0);

    // Store previous trial state and alpha for multi-phase approach
    auto prevTrialState = trialState;
    auto prevAlpha = alpha;

    // Main backtracking lambda
    auto backtrack = [&, this](bool isStrict) -> bool
    {
      std::string phaseLabel = isStrict ? "strict" : "relaxed";
      scalar_type zeta = isStrict ? one : zeta_;
      PRESSIOLOG_DEBUG("Start " + phaseLabel + " backtracking");

      while (true) {

        /*
         * Exit if alpha drops below 0.001; anything smaller
         * than this is probably unreasonable. Note that this
         * quantity doesn't depend on the dimension or magnitude
         * of the state.
         */
        if (std::abs(alpha) <= alphaLowerBound) {

          // If the strict phase fails, move on to relaxed phase
          if (isStrict) {
            trialState = prevTrialState;
            alpha = prevAlpha;
            PRESSIOLOG_INFO("Strict backtrack failed. Trying relaxed phase with alpha={}.", alpha);
            return false;
          }

          // If the relaxed phase fails, throw an error
          const auto msg = ": in BacktrackStrictlyDecreasingObjective: alpha = "
                         + std::to_string(alpha) + " <= " + std::to_string(alphaLowerBound)
                         + ", too small, exiting line search.";
          throw ::pressio::eh::LineSearchStepTooSmall(msg);
        }

        // During strict phase, store current values for (possible) future use
        if (isStrict) {
          prevAlpha = alpha;
          prevTrialState = trialState;
        }

        ::pressio::ops::update(trialState, zero, state, one, p_k, alpha);

        if (objective(trialState) <= zeta * objectiveValueAtCurrentNewtonStep){
          PRESSIOLOG_DEBUG(phaseLabel + " backtrack; condition satisfied with alpha= {:6e}", alpha);
          ::pressio::ops::update(state, one, p_k, alpha);
          return true;
        }

        // Conventional way to backtrack: alpha_l+1 = 0.5 * alpha_l
        alpha *= 0.5;
      }
    };

    auto strictBacktrack  = [&, this]() -> bool { return backtrack(true);  };
    auto relaxedBacktrack = [&, this]() -> bool { return backtrack(false); };

    // Try strict backtracking (zeta=1.0) first
    auto converged = strictBacktrack();
    if (!converged) (void)relaxedBacktrack();
  }
};

template<class RegistryType, class ObjF, class ScalarType, class StateType>
auto lm_gain_factor(RegistryType & reg,
		    ObjF & objective,
		    ScalarType objectiveValueAtCurrentNewtonStep,
		    StateType & cDiagH)
{
  using scalar_type = std::remove_const_t<ScalarType>;
  constexpr auto zero = static_cast<scalar_type>(0);
  constexpr auto one  = static_cast<scalar_type>(1);
  constexpr auto two  = static_cast<scalar_type>(2);

  const auto & state = reg.template get<StateTag>();
  const auto & correction  = reg.template get<CorrectionTag>();
  auto & trialState  = reg.template get<LineSearchTrialStateTag>();
  const auto & g = reg.template get<GradientTag>();
  const auto & H = reg.template get<LevenbergMarquardtUndampedHessianTag>();
  const auto & damp = reg.template get<LevenbergMarquardtDampingTag>();

  // numerator
  ::pressio::ops::update(trialState, zero, state, one, correction, one);
  const auto num = objective(trialState) - objectiveValueAtCurrentNewtonStep;

  // denominator
  const auto diagH = ::pressio::diagonal(H);
  ::pressio::ops::elementwise_multiply(one, correction, diagH, zero, cDiagH);
  const auto den1 = ::pressio::ops::dot(correction, cDiagH);
  const auto den2 = ::pressio::ops::dot(correction, g);
  const auto denom = (one/two)*(damp*den1 + den2);

  return num/denom;
}


template<class ScalarType, class StateType>
class LMSchedule1Updater
{
  using scalar_type = std::remove_const_t<ScalarType>;
  const scalar_type one = static_cast<scalar_type>(1);
  const scalar_type two = static_cast<scalar_type>(2);
  const scalar_type three = static_cast<scalar_type>(3);

  const scalar_type beta_     = two;
  const scalar_type gammaInv_ = one/three;
  const scalar_type p_ = three;
  const scalar_type tau_ = one;
  scalar_type nu_ = two;
  StateType cDiagH_;

public:
  LMSchedule1Updater(const StateType & stateIn) : cDiagH_(ops::clone(stateIn)){}

public:
  template<class RegType, class Objective>
  void operator()(RegType & reg,
		  Objective obj,
		  scalar_type objectiveValueAtCurrentNewtonStep)
  {
    PRESSIOLOG_DEBUG("nonlinsolver: lm1 update");

    auto & damp = reg.template get<LevenbergMarquardtDampingTag>();
    const auto & correction  = reg.template get<CorrectionTag>();
    auto & state = reg.template get<StateTag>();

    const scalar_type rho = lm_gain_factor(reg, obj, objectiveValueAtCurrentNewtonStep, cDiagH_);
    constexpr auto one  = static_cast<scalar_type>(1);
    constexpr auto two  = static_cast<scalar_type>(2);
    if (rho > 0){
      PRESSIOLOG_DEBUG("lm1 update: rho>0");
      ::pressio::ops::update(state, one, correction, one);
      const scalar_type a = one - (beta_ - one)*std::pow(two*rho - one, p_);
      damp *= std::max(gammaInv_, a);
      nu_ = beta_;
      PRESSIOLOG_DEBUG("lm1 update: rho>0, {}, {}", nu_, damp);
    }
    else{
      damp *= nu_;
      nu_ *= two;
      PRESSIOLOG_DEBUG("lm1 update: rho<=0, {}, {}", nu_, damp);
    }
  }
};

template<class ScalarType, class StateType>
class LMSchedule2Updater
{
  using scalar_type = std::remove_const_t<ScalarType>;

  const scalar_type rho1_	   = static_cast<scalar_type>(0.2);
  const scalar_type rho2_     = static_cast<scalar_type>(0.8);
  const scalar_type beta_	   = static_cast<scalar_type>(2);
  const scalar_type gammaInv_ = static_cast<scalar_type>(1)/static_cast<scalar_type>(3);
  const scalar_type tau_	   = static_cast<scalar_type>(1);
  StateType cDiagH_;

public:
  LMSchedule2Updater(const StateType & stateIn) : cDiagH_(ops::clone(stateIn)){}

public:
  template<class RegType, class Objective>
  void operator()(RegType & reg,
		  Objective obj,
		  scalar_type objectiveValueAtCurrentNewtonStep)
  {
    PRESSIOLOG_DEBUG("nonlinsolver: lm2 update");
    auto & damp = reg.template get<LevenbergMarquardtDampingTag>();
    const auto & correction  = reg.template get<CorrectionTag>();
    auto & state = reg.template get<StateTag>();

    constexpr auto one  = static_cast<scalar_type>(1);
    constexpr auto ten  = static_cast<scalar_type>(10);
    constexpr auto seven  = static_cast<scalar_type>(7);
    constexpr auto negSeven  = static_cast<scalar_type>(-1) * seven;
    const auto tenToSev  = std::pow(ten, seven);
    const auto tenToNegSev  = std::pow(ten, negSeven);

    const scalar_type rho = lm_gain_factor(reg, obj, objectiveValueAtCurrentNewtonStep, cDiagH_);
    if (rho < rho1_){
      damp = std::min(damp*beta_, tenToSev);
    }

    if (rho > rho2_){
      damp = std::max( tenToNegSev, damp*gammaInv_);
    }

    if (rho > 0){
      ::pressio::ops::update(state, one, correction, one);
    }
  }
};

}}}
#endif  // PRESSIOROM_SOLVERS_NONLINEAR_IMPL_UPDATERS_HPP_
