
#include <gtest/gtest.h>
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear_newton.hpp"
#include "./problems/problem2.hpp"

TEST(solvers_nonlinear, problem1A)
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  using namespace pressio;
  namespace pnonls = pressio::nonlinearsolvers;

  using problem_t  = solvers::test::Problem2;
  using state_t    = problem_t::state_type;
  using jacobian_t = problem_t::jacobian_type;

  problem_t sys;
  state_t y(2);
  y(0) = 0.3; y(1) = 0.4;

  using lin_solver_t = linearsolvers::Solver<linearsolvers::iterative::LSCG, jacobian_t>;
  lin_solver_t linearSolverObj;

  auto NonLinSolver = pnonls::create_newton_solver(sys, linearSolverObj);
  const auto updateMethod = pnonls::Update::BacktrackStrictlyDecreasingObjective;
  NonLinSolver.setUpdateCriterion(updateMethod);
  NonLinSolver.solve(y);

  const auto e1 = std::abs(y(0) - (1.93185165));
  const auto e2 = std::abs(y(1) - (0.51763809));
  ASSERT_TRUE(e1<1e-7);
  ASSERT_TRUE(e2<1e-7);
  std::cout << y << std::endl;

  PRESSIOLOG_FINALIZE();
}
