
#include <gtest/gtest.h>
#include "pressio/solvers.hpp"
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "testing_apps.hpp"

TEST(ode, implicit_bdf1_policy_default_created)
{
  using namespace pressio;
  using problem_t   = ode::testing::AppEigenB;
  using state_t = typename problem_t::state_type;
  problem_t problemObj;
  state_t y(problemObj.getInitCond());
  auto stepperObj = ode::create_implicit_stepper(ode::StepScheme::BDF1, problemObj);

  using jac_t = typename problem_t::jacobian_type;
  using lin_solver_t = linearsolvers::Solver<linearsolvers::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver = create_newton_solver(stepperObj,linSolverObj);

  ode::StepCount nSteps(2);
  double dt = 0.01;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, nSteps, NonLinSolver);
  std::cout << std::setprecision(14) << y << "\n";

  problemObj.analyticAdvanceBackEulerNSteps(dt, nSteps.get());
  EXPECT_DOUBLE_EQ(y(0), problemObj.y(0));
  EXPECT_DOUBLE_EQ(y(1), problemObj.y(1));
  EXPECT_DOUBLE_EQ(y(2), problemObj.y(2));
}

TEST(ode, implicit_bdf1_custom_policy)
{
  using namespace pressio;
  using problem_t = ode::testing::AppEigenB;
  using time_type = typename problem_t::independent_variable_type;
  using state_t = typename problem_t::state_type;
  using res_t = typename problem_t::rhs_type;
  using jac_t = typename problem_t::jacobian_type;

  problem_t problemObj;
  state_t y = problemObj.getInitCond();

  using rj_pol_t = ode::impl::ResidualJacobianStandardPolicy<
    problem_t, time_type, state_t, res_t, jac_t>;
  auto stepperObj = ode::create_implicit_stepper(ode::StepScheme::BDF1,
						 rj_pol_t(problemObj));

  using lin_solver_t = linearsolvers::Solver<linearsolvers::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;
  auto NonLinSolver = create_newton_solver(stepperObj,linSolverObj);

  ode::StepCount nSteps(2);
  double dt = 0.01;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, nSteps, NonLinSolver);
  std::cout << std::setprecision(14) << y << "\n";

  problemObj.analyticAdvanceBackEulerNSteps(dt, nSteps.get());
  EXPECT_NEAR(y(0), problemObj.y(0), 1e-15);
  EXPECT_NEAR(y(1), problemObj.y(1), 1e-15);
  EXPECT_NEAR(y(2), problemObj.y(2), 1e-15);
}

