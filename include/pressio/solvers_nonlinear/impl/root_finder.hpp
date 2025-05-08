
#ifndef SOLVERS_NONLINEAR_IMPL_ROOT_FINDER_HPP_
#define SOLVERS_NONLINEAR_IMPL_ROOT_FINDER_HPP_

#include <utility>

namespace pressio{
namespace nonlinearsolvers{
namespace impl{

template<
  class UserDefinedSystemType,
  class RegistryType,
  class ToleranceType,
  class NormDiagnosticsContainerType,
  class DiagnosticsLoggerType,
  class UpdaterType>
void root_solving_loop_impl(NewtonTag /*problemTag*/,
          const UserDefinedSystemType & system,
          RegistryType & reg,
          Stop stopEnumValue,
          ToleranceType stopTolerance,
          NormDiagnosticsContainerType & normDiagnostics,
          const DiagnosticsLoggerType & logger,
          int maxIters,
          UpdaterType && updater)
{

  using state_type = typename UserDefinedSystemType::state_type;

  auto objective = [&reg, &system](const state_type & stateIn){
    auto & r = reg.template get<ResidualTag>();
    system.residualAndJacobian(stateIn, r, {});
    return ::pressio::ops::norm2(r);
  };

  auto mustStop = [
		   &normDiag = std::as_const(normDiagnostics),
		   stopEnumValue, maxIters, stopTolerance]
    (int stepCount)
  {
    const Diagnostic stopDiag = stop_criterion_to_public_diagnostic(stopEnumValue);
    const bool hasConverged =
      is_absolute_diagnostic(stopDiag) ?
      normDiag[stopDiag].getAbsolute() < stopTolerance :
      normDiag[stopDiag].getRelative() < stopTolerance;
    const bool hasReachedMaxIters = stepCount >= maxIters;
    if (hasReachedMaxIters && !hasConverged) {
      writeNonlinearSolverTerminationFile(
        "Root Finder: Reached maximum number of iterations before convergence."
      );
    }
    switch (stopEnumValue){
    case Stop::AfterMaxIters:
      return hasReachedMaxIters;
    default:
      return hasConverged;
    };
  };

  int iStep = 0;
  while (++iStep <= maxIters){
    /* stage 1 */
    try{
      compute_residual_and_jacobian(reg, system);
    }
    catch (::pressio::eh::ResidualEvaluationFailureUnrecoverable const &e){
      PRESSIOLOG_ERROR(e.what());
      writeNonlinearSolverTerminationFile(e.what());
      throw ::pressio::eh::NonlinearSolveFailure();
    }
    catch (::pressio::eh::ResidualHasNans const &e){
      PRESSIOLOG_ERROR(e.what());
      writeNonlinearSolverTerminationFile(e.what());
      throw ::pressio::eh::NonlinearSolveFailure();
    }

    /* stage 2 */
    solve_newton_step(reg);

    /* stage 3 */
    std::for_each(normDiagnostics.begin(), normDiagnostics.end(),
      [&reg, iStep](auto & v){
        const bool isFirstIteration = iStep==1;
        compute_norm_internal_diagnostics(reg, isFirstIteration, v);
      });
    logger(iStep, normDiagnostics);

    /* stage 4*/
    if (mustStop(iStep)){
      PRESSIOLOG_DEBUG("nonlinsolver: stopping");
      break;
    }

    /* stage 5 */
    try{
      const auto currentObjValue =
	normDiagnostics[InternalDiagnostic::residualAbsoluteRelativel2Norm].getAbsolute();
        updater(reg, objective, currentObjValue);
    }
    catch (::pressio::eh::LineSearchStepTooSmall const &e) {
      // nicely exit the solve
      PRESSIOLOG_WARNING(e.what());
      writeNonlinearSolverTerminationFile(e.what());
      break;
    }
    catch (::pressio::eh::LineSearchObjFunctionChangeTooSmall const &e) {
      // nicely exit the solve
      PRESSIOLOG_WARNING(e.what());
      writeNonlinearSolverTerminationFile(e.what());
      break;
    }
  }
}

  /*
    this wrapper is needed because the user-defined problem needs
    to know what is the current nonlinear state to use for computing
    the action of the jacobian
  */
template<class UserSystemT>
class Wrapper
{
  using state_type = typename UserSystemT::state_type;
  UserSystemT const * m_system = nullptr;
  state_type const * m_currentNonLinearState = nullptr;
  mutable state_type m_auxiliary;
  int m_numJacRows = {};
  int m_numJacCols = {};

public:
  // needed this by the gmres solver
  using scalar_type = scalar_trait_t<state_type>;

  Wrapper(UserSystemT const * system, int numRows, int numCols)
    : m_system(system),
      m_auxiliary(system->createState()),
      m_numJacRows(numRows), m_numJacCols(numCols)
  {}

  int rows() const { return m_numJacRows; }
  int cols() const { return m_numJacCols; }

  void setCurrentNonLinearState(typename UserSystemT::state_type const & newState){
    m_currentNonLinearState = &newState;
  }

  template<class OperandT, class ResultT>
  void applyAndAddTo(OperandT const & operand, ResultT & out) const {
    // we need to compute out = out + jacobian * operand
    // since user computes jacobian * operand, we need to do one extra step here
    // so we do:
    // m_auxiliary = jacobian * operand
    // out += m_auxiliary

    pressio::ops::set_zero(m_auxiliary);
    m_system->applyJacobian(*m_currentNonLinearState, operand, m_auxiliary);
    pressio::ops::update(out, 1, m_auxiliary, 1);
  }
};

template<
  class UserDefinedSystemType,
  class RegistryType,
  class ToleranceType,
  class NormDiagnosticsContainerType,
  class DiagnosticsLoggerType,
  class UpdaterType>
void root_solving_loop_impl(MatrixFreeNewtonTag /*problemTag*/,
          const UserDefinedSystemType & system,
          RegistryType & reg,
          Stop stopEnumValue,
          ToleranceType stopTolerance,
          NormDiagnosticsContainerType & normDiagnostics,
          const DiagnosticsLoggerType & logger,
          int maxIters,
          UpdaterType && updater)
{

  using state_type = typename UserDefinedSystemType::state_type;

  auto objective = [&reg, &system](const state_type & stateIn){
    auto & r = reg.template get<ResidualTag>();
    system.residual(stateIn, r);
    return ::pressio::ops::norm2(r);
  };

  auto mustStop = [
		   &normDiag = std::as_const(normDiagnostics),
		   stopEnumValue, maxIters, stopTolerance]
    (int stepCount)
  {
    const Diagnostic stopDiag = stop_criterion_to_public_diagnostic(stopEnumValue);
    const bool hasConverged =
      is_absolute_diagnostic(stopDiag) ?
      normDiag[stopDiag].getAbsolute() < stopTolerance :
      normDiag[stopDiag].getRelative() < stopTolerance;
    const bool hasReachedMaxIters = stepCount >= maxIters;
    if (hasReachedMaxIters && !hasConverged) {
      writeNonlinearSolverTerminationFile(
        "Root Finder: Reached maximum number of iterations before convergence."
      );
    }
    switch (stopEnumValue){
    case Stop::AfterMaxIters:
      return hasReachedMaxIters;
    default:
      return hasConverged;
    };
  };

  using wrapper_t = Wrapper<UserDefinedSystemType>;
  using linear_solver_t =  pressio::linearsolvers::Solver<
    typename RegistryType::extended_registry::linear_solver_tag,
    wrapper_t>;
  linear_solver_t linearSolver;

  const auto & currentState = reg.template get<StateTag>();
  const auto & r = reg.template get<ResidualTag>();
  wrapper_t wrapper(&system,
		    pressio::ops::extent(r, 0),
		    pressio::ops::extent(currentState, 0));

  int iStep = 0;
  while (++iStep <= maxIters)
  {
    /* stage 1 */
    try{
      compute_residual(reg, system);
    }
    catch (::pressio::eh::ResidualEvaluationFailureUnrecoverable const &e){
      PRESSIOLOG_ERROR(e.what());
      writeNonlinearSolverTerminationFile(e.what());
      throw ::pressio::eh::NonlinearSolveFailure();
    }
    catch (::pressio::eh::ResidualHasNans const &e){
      PRESSIOLOG_ERROR(e.what());
      writeNonlinearSolverTerminationFile(e.what());
      throw ::pressio::eh::NonlinearSolveFailure();
    }

    /* stage 2, matrix-free newton step */
    wrapper.setCurrentNonLinearState(currentState);
    auto & c = reg.template get<CorrectionTag>();
    linearSolver.solve(wrapper, r, c);

    // scale by -1 for sign convention
    using c_t = mpl::remove_cvref_t<decltype(c)>;
    using scalar_type = typename ::pressio::Traits<c_t>::scalar_type;
    pressio::ops::scale(c, static_cast<scalar_type>(-1));

    /* stage 3 */
    std::for_each(normDiagnostics.begin(), normDiagnostics.end(),
      [&reg, iStep](auto & v){
        [[maybe_unused]] const bool isFirstIteration = iStep==1;
	compute_norm_internal_diagnostics(reg, isFirstIteration, v);
      });
    logger(iStep, normDiagnostics);

    /* stage 4*/
    if (mustStop(iStep)){
      PRESSIOLOG_DEBUG("nonlinsolver: stopping");
      break;
    }

    /* stage 5 */
    try{
      const auto currentObjValue =
	normDiagnostics[InternalDiagnostic::residualAbsoluteRelativel2Norm].getAbsolute();
      updater(reg, objective, currentObjValue);
    }
    catch (::pressio::eh::LineSearchStepTooSmall const &e) {
      // nicely exit the solve
      PRESSIOLOG_WARNING(e.what());
      writeNonlinearSolverTerminationFile(e.what());
      break;
    }
    catch (::pressio::eh::LineSearchObjFunctionChangeTooSmall const &e) {
      // nicely exit the solve
      PRESSIOLOG_WARNING(e.what());
      writeNonlinearSolverTerminationFile(e.what());
      break;
    }
  }
}


template<class Tag, class StateType, class RegistryType, class NormValueType>
class RootFinder : public RegistryType
{
  Tag tag_;
  int maxIters_ = 100;
  Stop stopEnValue_ = Stop::WhenAbsolutel2NormOfCorrectionBelowTolerance;
  NormValueType stopTolerance_ = 0.000001;
  Update updateEnValue_ = Update::Standard;

  using norm_diagnostics_container = DiagnosticsContainer<
    InternalDiagnosticDataWithAbsoluteRelativeTracking<NormValueType> >;
  norm_diagnostics_container normDiagnostics_;
  DiagnosticsLogger diagnosticsLogger_ = {};
  std::optional<std::vector<scalar_trait_t<StateType> > > parameters_;

public:
  template<class ...Args>
  RootFinder(Tag tagIn,
	     const std::vector<Diagnostic> & diags,
	     Args && ...args)
    : RegistryType(std::forward<Args>(args)...),
      tag_(tagIn),
      normDiagnostics_(diags)
  {

    // currently we don't have the diagonostics stuff all flushed out
    // so we limit it to work for a specific case
    const auto & publicDiags = normDiagnostics_.publicNames();
    assert(publicDiags.size() == 4);
    assert(publicDiags[0] == Diagnostic::residualAbsolutel2Norm);
    assert(publicDiags[1] == Diagnostic::residualRelativel2Norm);
    assert(publicDiags[2] == Diagnostic::correctionAbsolutel2Norm);
    assert(publicDiags[3] == Diagnostic::correctionRelativel2Norm);

    // check the stop criterion uses a metric already supported in the diagonstics
    const auto stopMetric = stop_criterion_to_public_diagnostic(stopEnValue_);
    normDiagnostics_.addIfUnsupported(stopMetric);
    // need to reset the logger since the names might have changed
    diagnosticsLogger_.resetFor(publicDiags);
  }

  // query/set update criterion
  Update currentUpdateCriterion() const   { return updateEnValue_; }
  void setUpdateCriterion(Update value) { updateEnValue_ = value; }

  // query/set stop criterion, tolerance, line search params
  Stop currentStopCriterion() const          { return stopEnValue_; }
  void setStopCriterion(Stop value)	     { stopEnValue_ = value; }
  void setStopTolerance(NormValueType value) { stopTolerance_ = value; }
  void setMaxIterations(int newMax)          { maxIters_ = newMax; }
  auto & getLineSearchParameters(){
    return parameters_;
  }
  auto addLineSearchParameter(scalar_trait_t<StateType> new_param) {
    if (!parameters_) {
      parameters_ = std::vector<scalar_trait_t<StateType>>{};
    }
    parameters_->push_back(new_param);
  }

  template<class SystemType, class _Tag = Tag>
  std::enable_if_t< std::is_same<_Tag, NewtonTag>::value >
  solve(const SystemType & system, StateType & solutionInOut)
  {
    if (updateEnValue_ == Update::Standard)
    {
      auto extReg = reference_capture_registry_and_extend_with<
	StateTag, StateType &>(*this, solutionInOut);
      root_solving_loop_impl(tag_, system, extReg, stopEnValue_, stopTolerance_,
			     normDiagnostics_, diagnosticsLogger_, maxIters_,
			     DefaultUpdater());

    }
    else if (updateEnValue_ == Update::BacktrackStrictlyDecreasingObjective)
    {
      auto extReg = reference_capture_registry_and_extend_with<
	StateTag, LineSearchTrialStateTag,
	StateType &, StateType>(*this, solutionInOut, system.createState());

      root_solving_loop_impl(tag_, system, extReg, stopEnValue_, stopTolerance_,
			     normDiagnostics_, diagnosticsLogger_, maxIters_,
			     BacktrackStrictlyDecreasingObjectiveUpdater<scalar_trait_t<StateType>>(parameters_));
    }
    else{
      throw std::runtime_error("Invalid criterion");
    }
  }

  template<class SystemType, class _Tag = Tag>
  std::enable_if_t< std::is_same<_Tag, MatrixFreeNewtonTag>::value >
  solve(const SystemType & system, StateType & solutionInOut)
  {
    assert(updateEnValue_ == Update::Standard);

    auto extReg = reference_capture_registry_and_extend_with<
      StateTag, StateType &>(*this, solutionInOut);
    root_solving_loop_impl(tag_, system, extReg, stopEnValue_, stopTolerance_,
			   normDiagnostics_, diagnosticsLogger_, maxIters_,
			   DefaultUpdater());
  }

  // this method can be used when the solver is applied
  // to the same system used for constructing it
  void solve(StateType & solutionInOut)
  {
    auto * system = this->template get<SystemTag>();
    assert(system != nullptr);
    this->solve(*system, solutionInOut);
  }
};


}}}
#endif
