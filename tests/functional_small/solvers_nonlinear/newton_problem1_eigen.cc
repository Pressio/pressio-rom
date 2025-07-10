
#include <gtest/gtest.h>
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear_newton.hpp"
#include "./problems/problem1.hpp"

template<class SystemType>
void run_impl(int reps, bool logOn = false, bool callSolveWithJustState = true)
{
  if (logOn){
    PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);
  }

  using namespace pressio;
  using problem_t  = SystemType;
  using state_t    = typename problem_t::state_type;
  using jacobian_t = typename problem_t::jacobian_type;

  using lin_solver_t = linearsolvers::Solver<linearsolvers::iterative::LSCG, jacobian_t>;
  lin_solver_t linearSolverObj;

  problem_t sys;
  state_t y(2);
  auto nonLinSolver = create_newton_solver(sys, linearSolverObj);

  if (reps){
    y(0) = 0.001; y(1) = 0.0001;

    if (callSolveWithJustState){
      nonLinSolver.solve(y);
    }else{
      nonLinSolver.solve(sys, y);
    }
    const auto e1 = std::abs(y(0) - (1.));
    const auto e2 = std::abs(y(1) - (0.));
    std::cout << y << std::endl;
    ASSERT_TRUE(e1<1e-8);
    ASSERT_TRUE(e2<1e-8);
  }
  else{

    for (int i=0; i<reps; ++i){
      y(0) = 0.001; y(1) = 0.0001;
      if (callSolveWithJustState){
	nonLinSolver.solve(y);
      }else{
	nonLinSolver.solve(sys, y);
      }

      const auto e1 = std::abs(y(0) - (1.));
      const auto e2 = std::abs(y(1) - (0.));
      ASSERT_TRUE(e1<1e-8);
      ASSERT_TRUE(e2<1e-8);
    }
  }

  if (logOn){
    PRESSIOLOG_FINALIZE();
  }
}

TEST(solvers_nonlinear, problem1){
  run_impl<pressio::solvers::test::Problem1>(1, true, false);
}

TEST(solvers_nonlinear, problem1_repeated_solve){
  run_impl<pressio::solvers::test::Problem1>(100, false, false);
}

TEST(solvers_nonlinear, problem1_call_solve_with_only_state){
  run_impl<pressio::solvers::test::Problem1>(1, true, true);
}

TEST(solvers_nonlinear, problem1_repeated_solve_call_solve_with_only_state){
  run_impl<pressio::solvers::test::Problem1>(100, false, true);
}

TEST(solvers_nonlinear, move)
{
  using namespace pressio;
  using problem_t  = pressio::solvers::test::Problem1;
  using state_t    = typename problem_t::state_type;
  using jacobian_t = typename problem_t::jacobian_type;

  using lin_solver_t = linearsolvers::Solver<linearsolvers::iterative::LSCG, jacobian_t>;
  lin_solver_t linearSolverObj;

  problem_t sys;
  auto dosolve = [&](state_t & y, auto & nls){
    y(0) = 0.001; y(1) = 0.0001;
    nls.solve(sys, y);
    const auto e1 = std::abs(y(0) - (1.));
    const auto e2 = std::abs(y(1) - (0.));
    std::cout << y << std::endl;
    ASSERT_TRUE(e1<1e-8);
    ASSERT_TRUE(e2<1e-8);
  };

  state_t y(2);
  auto nls = create_newton_solver(sys, linearSolverObj);
  dosolve(y, nls);

  auto nls2 = std::move(nls);
  dosolve(y, nls2);
}
