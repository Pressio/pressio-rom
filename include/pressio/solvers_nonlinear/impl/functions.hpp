
#ifndef PRESSIO_SOLVERS_NONLINEAR_IMPL_FUNCTIONS_HPP_
#define PRESSIO_SOLVERS_NONLINEAR_IMPL_FUNCTIONS_HPP_

#include <fstream>

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

#ifdef PRESSIO_ENABLE_CXX20
  template<class RegistryType, class StateType, class SystemType>
  requires NonlinearSystemFusingResidualAndJacobian<SystemType>
#else
  template<
    class RegistryType, class StateType, class SystemType,
    std::enable_if_t<NonlinearSystemFusingResidualAndJacobian<SystemType>::value, int> = 0
  >
#endif
void compute_residual(RegistryType & reg,
		      const StateType & state,
		      const SystemType & system)
{
  auto & r = reg.template get<ResidualTag>();
  system.residualAndJacobian(state, r, {});
}

#ifdef PRESSIO_ENABLE_CXX20
template<class RegistryType, class StateType, class SystemType>
requires NonlinearSystem<SystemType>
#else
template<
  class RegistryType, class StateType, class SystemType,
  std::enable_if_t< NonlinearSystem<SystemType>::value, int> = 0
  >
#endif
void compute_residual(RegistryType & reg,
		      const StateType & state,
		      const SystemType & system)
{
  auto & r = reg.template get<ResidualTag>();
  system.residual(state, r);
}

#ifdef PRESSIO_ENABLE_CXX20
template<class RegistryType, class SystemType>
requires NonlinearSystem<SystemType>
#else
template<
  class RegistryType, class SystemType,
  std::enable_if_t< NonlinearSystem<SystemType>::value, int> = 0
  >
#endif
void compute_residual(RegistryType & reg,
		      const SystemType & system)
{
  const auto & state = reg.template get<StateTag>();
  auto & r = reg.template get<ResidualTag>();
  system.residual(state, r);
}

template<class RegistryType, class SystemType>
void compute_residual_and_jacobian(RegistryType & reg,
				   const SystemType & system)
{
  const auto & state = reg.template get<StateTag>();
  auto & r = reg.template get<ResidualTag>();
  auto & j = reg.template get<JacobianTag>();
  using j_t = typename SystemType::jacobian_type;
  system.residualAndJacobian(state, r, std::optional<j_t*>{&j});
}

template<class RegistryType>
void compute_gradient(RegistryType & reg)
{
  constexpr auto pT  = ::pressio::transpose();
  const auto & r = reg.template get<ResidualTag>();
  const auto & j = reg.template get<JacobianTag>();
  auto & g = reg.template get<GradientTag>();
  // g = J_r^T r
  ::pressio::ops::product(pT, 1, j, r, 0, g);
}

template<class T>
auto compute_half_sum_of_squares(const T & operand)
{
  static_assert(Traits<T>::rank == 1, "");
  const auto normVal = ::pressio::ops::norm2(operand);
  using sc_type = mpl::remove_cvref_t<decltype(normVal)>;
  constexpr auto one  = static_cast<sc_type>(1);
  constexpr auto two  = static_cast<sc_type>(2);
  return std::pow(normVal, two)*(one/two);
}

template<class RegistryType, class StateType, class SystemType>
auto compute_nonlinearls_objective(GaussNewtonNormalEqTag /*tag*/,
				   RegistryType & reg,
				   const StateType & state,
				   const SystemType & system)
{
  compute_residual(reg, state, system);
  auto & r = reg.template get<ResidualTag>();
  return compute_half_sum_of_squares(r);
}

template<class RegistryType, class StateType, class SystemType>
auto compute_nonlinearls_objective(GaussNewtonQrTag /*tag*/,
				   RegistryType & reg,
				   const StateType & state,
				   const SystemType & system)
{
  compute_residual(reg, state, system);
  auto & r = reg.template get<ResidualTag>();
  return compute_half_sum_of_squares(r);
}

template<class RegistryType, class StateType, class SystemType>
auto compute_nonlinearls_objective(LevenbergMarquardtNormalEqTag /*tag*/,
				   RegistryType & reg,
				   const StateType & state,
				   const SystemType & system)
{
  compute_residual(reg, state, system);
  auto & r = reg.template get<ResidualTag>();
  return compute_half_sum_of_squares(r);
}

template<class RegistryType, class StateType, class SystemType>
auto compute_nonlinearls_objective(WeightedGaussNewtonNormalEqTag /*tag*/,
				   RegistryType & reg,
				   const StateType & state,
				   const SystemType & system)
{
  const auto & W = reg.template get<WeightingOperatorTag>();
  auto & r  = reg.template get<ResidualTag>();
  auto & Wr = reg.template get<WeightedResidualTag>();
  compute_residual(reg, state, system);
  (*W)(r, Wr);

  const auto v = ::pressio::ops::dot(r, Wr);
  using sc_t = mpl::remove_cvref_t< decltype(v) >;
  return v * (static_cast<sc_t>(1) / static_cast<sc_t>(2));
}

template<class RegistryType, class StateType, class SystemType>
auto compute_nonlinearls_objective(CompactWeightedGaussNewtonNormalEqTag /*tag*/,
				   RegistryType & reg,
				   const StateType & state,
				   const SystemType & system)
{
  const auto & W = reg.template get<WeightingOperatorTag>();
  auto & r  = reg.template get<ResidualTag>();
  auto & Wr = reg.template get<WeightedResidualTag>();
  compute_residual(reg, state, system);
  (*W)(r, Wr);

  const auto v = ::pressio::ops::dot(Wr, Wr);
  using sc_t = mpl::remove_cvref_t< decltype(v) >;
  return v * (static_cast<sc_t>(1) / static_cast<sc_t>(2));
}

#ifdef PRESSIO_ENABLE_CXX20
template<class RegistryType, class SystemType>
requires RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>
#else
template<
  class RegistryType, class SystemType,
  std::enable_if_t<
    RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>::value,
    int> = 0
  >
#endif
auto compute_nonlinearls_operators_and_objective(GaussNewtonNormalEqTag /*tag*/,
						 RegistryType & reg,
						 const SystemType & system)
{

  compute_residual_and_jacobian(reg, system);

  const auto & r = reg.template get<ResidualTag>();
  const auto & J = reg.template get<JacobianTag>();
  auto & g = reg.template get<GradientTag>();
  auto & H = reg.template get<HessianTag>();
  constexpr auto pT  = ::pressio::transpose();
  constexpr auto pnT = ::pressio::nontranspose();
  // H = J_r^T J_r, g = J_r^T r
  ::pressio::ops::product(pT, pnT, 1, J, 0, H);
  ::pressio::ops::product(pT, 1, J, r, 0, g);

  return compute_half_sum_of_squares(r);
}

#ifdef PRESSIO_ENABLE_CXX20
template<class RegistryType, class SystemType>
requires RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>
#else
template<
  class RegistryType, class SystemType,
  std::enable_if_t<
    RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>::value,
    int> = 0
  >
#endif
auto compute_nonlinearls_operators_and_objective(GaussNewtonQrTag /*tag*/,
						 RegistryType & reg,
						 const SystemType & system)
{

  compute_residual_and_jacobian(reg, system);

  const auto & r = reg.template get<ResidualTag>();
  const auto & J = reg.template get<JacobianTag>();
  auto & g = reg.template get<GradientTag>();
  constexpr auto pT  = ::pressio::transpose();
  // g = J_r^T r
  ::pressio::ops::product(pT, 1, J, r, 0, g);

  return compute_half_sum_of_squares(r);
}

#ifdef PRESSIO_ENABLE_CXX20
template<class RegistryType, class SystemType>
requires RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>
#else
template<
  class RegistryType, class SystemType,
  std::enable_if_t<
    RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>::value,
    int> = 0
  >
#endif
auto compute_nonlinearls_operators_and_objective(WeightedGaussNewtonNormalEqTag /*tag*/,
						 RegistryType & reg,
						 const SystemType & system)
{
  compute_residual_and_jacobian(reg, system);

  constexpr auto pT  = ::pressio::transpose();
  constexpr auto pnT = ::pressio::nontranspose();
  const auto & W = reg.template get<WeightingOperatorTag>();
  const auto & r = reg.template get<ResidualTag>();
  const auto & J = reg.template get<JacobianTag>();
  auto & Wr = reg.template get<WeightedResidualTag>();
  auto & WJ = reg.template get<WeightedJacobianTag>();
  auto & g  = reg.template get<GradientTag>();
  auto & H  = reg.template get<HessianTag>();

  (*W)(r, Wr);
  (*W)(J, WJ);
  ::pressio::ops::product(pT, pnT, 1, J, WJ, 0, H);
  ::pressio::ops::product(pT, 1, J, Wr, 0, g);

  using sc_t = scalar_trait_t<typename SystemType::state_type>;
  const auto v = ::pressio::ops::dot(r, Wr);
  return v * (static_cast<sc_t>(1) / static_cast<sc_t>(2));
}

/* Special case of weighted Gauss Newton, just changing the action of W such that
    H = (J^T_r * W^T) * (W * J_r)
    g = (J^T_r * W^T) * (W * r)
  In instances where W is shape [M x N] and M << N, this results in a significant memory usage reduction
  Particularly useful for GNAT on large sample meshes,
    where M is the number of modes and N is the number of sampling points
*/
#ifdef PRESSIO_ENABLE_CXX20
template<class RegistryType, class SystemType>
requires RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>
#else
template<
  class RegistryType, class SystemType,
  std::enable_if_t<
    RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>::value,
    int> = 0
  >
#endif
auto compute_nonlinearls_operators_and_objective(CompactWeightedGaussNewtonNormalEqTag /*tag*/,
						 RegistryType & reg,
						 const SystemType & system)
{
  compute_residual_and_jacobian(reg, system);

  constexpr auto pT  = ::pressio::transpose();
  constexpr auto pnT = ::pressio::nontranspose();
  const auto & W = reg.template get<WeightingOperatorTag>();
  const auto & r = reg.template get<ResidualTag>();
  const auto & J = reg.template get<JacobianTag>();
  auto & Wr = reg.template get<WeightedResidualTag>();
  auto & WJ = reg.template get<WeightedJacobianTag>();
  auto & g  = reg.template get<GradientTag>();
  auto & H  = reg.template get<HessianTag>();

  (*W)(r, Wr);
  (*W)(J, WJ);
  ::pressio::ops::product(pT, pnT, 1, WJ, WJ, 0, H);
  ::pressio::ops::product(pT, 1, WJ, Wr, 0, g);

  using sc_t = scalar_trait_t<typename SystemType::state_type>;
  const auto v = ::pressio::ops::dot(Wr, Wr);
  return v * (static_cast<sc_t>(1) / static_cast<sc_t>(2));
}


#ifdef PRESSIO_ENABLE_CXX20
template<class RegistryType, class SystemType>
requires RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>
#else
template<
  class RegistryType, class SystemType,
  std::enable_if_t<
    RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>::value,
    int> = 0
  >
#endif
auto compute_nonlinearls_operators_and_objective(LevenbergMarquardtNormalEqTag /*tag*/,
						 RegistryType & reg,
						 const SystemType & system)
{
  compute_residual_and_jacobian(reg, system);

  constexpr auto pT  = ::pressio::transpose();
  constexpr auto pnT = ::pressio::nontranspose();
  const auto & r = reg.template get<ResidualTag>();
  const auto & J = reg.template get<JacobianTag>();
  auto & g = reg.template get<GradientTag>();
  auto & H = reg.template get<LevenbergMarquardtUndampedHessianTag>();
  auto & scaledH = reg.template get<HessianTag>();
  const auto & damp = reg.template get<LevenbergMarquardtDampingTag>();

  ::pressio::ops::product(pT, pnT, 1, J, 0, H);
  ::pressio::ops::product(pT, 1, J, r, 0, g);

  // compute scaledH = H + mu*diagonal(H)
  ::pressio::ops::deep_copy(scaledH, H);
  const auto diagH = ::pressio::diagonal(H);
  auto diaglmH = ::pressio::diagonal(scaledH);
  ::pressio::ops::update(diaglmH, 1, diagH, damp);

  return compute_half_sum_of_squares(r);
}

template<class RegistryType>
void solve_newton_step(RegistryType & reg)
{
  /* Newton correction solves: J_r delta = - r
     where: delta = x_k+1 - x_k
     so we solve: J_r (-delta) = r and then scale by -1
  */

  const auto & r = reg.template get<ResidualTag>();
  const auto & J = reg.template get<JacobianTag>();
  auto & c = reg.template get<CorrectionTag>();
  auto & solver = reg.template get<InnerSolverTag>();
  // solve J_r correction = r
  solver->solve(J, r, c);
  // scale by -1 for sign convention
  using c_t = mpl::remove_cvref_t<decltype(c)>;
  using scalar_type = typename ::pressio::Traits<c_t>::scalar_type;
  pressio::ops::scale(c, static_cast<scalar_type>(-1));
}

template<class RegistryType>
void solve_hessian_gradient_linear_system(RegistryType & reg)
{
  const auto & g = reg.template get<GradientTag>();
  const auto & H = reg.template get<HessianTag>();
  auto & c = reg.template get<CorrectionTag>();
  auto & solver = reg.template get<InnerSolverTag>();
  solver->solve(H, g, c);
}

template<class RegistryType>
void compute_correction(GaussNewtonNormalEqTag /*tag*/,
			RegistryType & reg)
{
  /* Gauss-newton with normal eq, we are solving:
       J_r^T*J_r  (x_k+1 - x_k) = - J_r^T*r
     which we can write as:  H c = g
       H = J_r^T*J_r
       g = J_r^T r
       c = x_k - x_k+1

     IMPORTANT: since we define the correction as:
       c = x_k+1 - x_k,
     we need to rescale c by -1 after the solve
  */
  solve_hessian_gradient_linear_system(reg);
  auto & c = reg.template get<CorrectionTag>();
  ::pressio::ops::scale(c, -1);
}

template<class RegistryType>
void compute_correction(WeightedGaussNewtonNormalEqTag /*tag*/,
			RegistryType & reg)
{
  // this is same as regular GN since we solve H delta = g
  solve_hessian_gradient_linear_system(reg);
  auto & c = reg.template get<CorrectionTag>();
  ::pressio::ops::scale(c, -1);
}

template<class RegistryType>
void compute_correction(CompactWeightedGaussNewtonNormalEqTag /*tag*/,
			RegistryType & reg)
{
  // this is same as regular GN since we solve H delta = g
  solve_hessian_gradient_linear_system(reg);
  auto & c = reg.template get<CorrectionTag>();
  ::pressio::ops::scale(c, -1);
}

template<class RegistryType>
void compute_correction(LevenbergMarquardtNormalEqTag /*tag*/,
			RegistryType & reg)
{
  /*
    For LM we are solving: H c = g
    where:
       H = J_r^T*J_r + lambda*diagonal(J_r^T J_r)
       g = J_r^T r
       c = x_k - x_k+1

    IMPORTANT: since we define the correction as:
       c = x_k+1 - x_k,
    we need to rescale c by -1 after
  */
  solve_hessian_gradient_linear_system(reg);
  auto & c = reg.template get<CorrectionTag>();
  ::pressio::ops::scale(c, -1);
}

template<class RegistryType>
void compute_correction(GaussNewtonQrTag /*tag*/,
			RegistryType & reg)
{
  /*
    see: https://en.wikipedia.org/wiki/Non-linear_least_squares#QR_decomposition
    but careful that here we use J to refer to jacobian wrt r,
    which is the negative of J_f

     IMPORTANT: since we define the correction as:
       c = x_k+1 - x_k,
     we need to rescale c by -1 after the solve
  */

  const auto & r  = reg.template get<ResidualTag>();
  const auto & J = reg.template get<JacobianTag>();
  auto & c = reg.template get<CorrectionTag>();
  auto & QTr = reg.template get<QTransposeResidualTag>();
  auto & solver = reg.template get<InnerSolverTag>();

  // factorize J = QR
  solver->computeThin(J);

  // compute Q^T r
  solver->applyQTranspose(r, QTr);

  // solve Rfactor c = Q^T r
  solver->solve(QTr, c);

  // rescale as said above
  ::pressio::ops::scale(c, -1);
}

// =====================================================

template<class Tag, class RegistryType>
void reset_for_new_solve_loop(Tag /*tag*/, RegistryType & reg){
  // no op
}

template<class RegistryType>
void reset_for_new_solve_loop(LevenbergMarquardtNormalEqTag /*tagJ*/,
			      RegistryType & reg)
{
  auto & damp = reg.template get<LevenbergMarquardtDampingTag>();
  damp = 1;
}

// =====================================================

template<class Tag, class T, class RegistryType>
std::enable_if_t< !RegistryType::template contains<Tag>() >
compute_norm2_if_tag_if_present(Tag /*t*/,
				const RegistryType & reg,
				bool isInitial,
				InternalDiagnosticDataWithAbsoluteRelativeTracking<T> & metric)
{ /* noop */ }

template<class Tag, class T, class RegistryType>
std::enable_if_t< RegistryType::template contains<Tag>() >
compute_norm2_if_tag_if_present(Tag /*t*/,
				const RegistryType & reg,
				bool isInitial,
				InternalDiagnosticDataWithAbsoluteRelativeTracking<T> & metric)
{
  const auto & operand = reg.template get<Tag>();
  const auto value = ops::norm2(operand);
  metric.update(value, isInitial);
}

template<class T, class RegistryType>
void compute_norm_internal_diagnostics(const RegistryType & reg,
				       bool isInitial,
				       InternalDiagnosticDataWithAbsoluteRelativeTracking<T> & metric)
{

  switch(metric.name())
  {
    case InternalDiagnostic::residualAbsoluteRelativel2Norm:
      compute_norm2_if_tag_if_present(ResidualTag{}, reg, isInitial, metric);
      break;

    case InternalDiagnostic::correctionAbsoluteRelativel2Norm:
      compute_norm2_if_tag_if_present(CorrectionTag{}, reg, isInitial, metric);
      break;

    case InternalDiagnostic::gradientAbsoluteRelativel2Norm:
      compute_norm2_if_tag_if_present(GradientTag{}, reg, isInitial, metric);
      break;

    default: return;
  };//end switch
}

// =====================================================

// Write a text file with the error that caused the
// nonlinear solver to exit without converging.
inline void writeNonlinearSolverTerminationFile(const std::string terminationString) {
  int rank = 0;
#if defined PRESSIO_ENABLE_TPL_MPI
  int flag = 0; MPI_Initialized( &flag );
  if (flag == 1) MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if (rank == 0) {
    const std::string fileName = "nonlinear_solver_failed.txt";
    std::ofstream file; file.open(fileName);
    file << terminationString << std::endl;
    file.close();
  }
}

}}}
#endif  // PRESSIO_SOLVERS_NONLINEAR_IMPL_FUNCTIONS_HPP_
