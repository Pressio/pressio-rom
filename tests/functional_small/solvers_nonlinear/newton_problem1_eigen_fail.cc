
#include <gtest/gtest.h>
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear_newton.hpp"
#include "./problems/problem1.hpp"

#include <filesystem>

/**
 * Tests all failure conditions RootFinder nonlinear solver.
 * Any failure to converge should write out "nonlinear_solver_fail.txt"
 * with the corresponding error message.
 */

enum class FailType {
        ResidualEvaluationFailureUnrecoverable,
        ResidualHasNans,
        MaximumIterations,
        LineSearchStepTooSmall,
        LineSearchObjFunctionChangeTooSmall
};

const std::unordered_map<FailType, std::string> failMsgs = {
    {
        FailType::ResidualEvaluationFailureUnrecoverable,
        std::string(pressio::eh::ResidualEvaluationFailureUnrecoverable().what())
    },
    {
        FailType::ResidualHasNans,
        std::string(pressio::eh::ResidualHasNans().what())
    },
    {
        FailType::MaximumIterations,
        "Root Finder: Reached maximum number of iterations before convergence."
    },
    {
        FailType::LineSearchStepTooSmall,
        std::string(pressio::eh::LineSearchStepTooSmall().what())
    },
    {
        FailType::LineSearchObjFunctionChangeTooSmall,
        std::string(pressio::eh::LineSearchObjFunctionChangeTooSmall().what())
    }
};

template <FailType failType>
struct AdaptedProblem1 : pressio::solvers::test::Problem1
{
  void residualAndJacobian(
    const state_type& x,
		residual_type& res,
		std::optional<jacobian_type*> Jin) const
  {
    if constexpr (failType == FailType::ResidualHasNans) {
      // This exception is never explicitly thrown in Pressio,
      // so we have to throw it manually here.
      throw pressio::eh::ResidualHasNans();
    } else if constexpr (failType == FailType::ResidualEvaluationFailureUnrecoverable) {
      // This exception is only thrown if it catches either of
      // the following:
      //    VelocityFailureUnrecoverable
      //    DiscreteTimeResidualFailureUnrecoverable
      // However, those are never explicitly thrown, and so
      // we have to manually throw the exception here.
      throw pressio::eh::ResidualEvaluationFailureUnrecoverable();
    } else if constexpr (failType == FailType::LineSearchObjFunctionChangeTooSmall) {
      // Makes the difference between new and old objective values
      // negligible, which should trigger the exception
      res(0) = 1e-14 * x(0);
      res(1) = 1e-14 * x(1);
    } else {
      // Standard values from Problem1 test.
      res(0) =  x(0)*x(0)*x(0) + x(1) - 1.0;
      res(1) = -x(0) + x(1)*x(1)*x(1) + 1.0;
    }

    if (Jin){
      auto & jac = *Jin.value();
      jac.coeffRef(0, 0) = 3.0*x(0)*x(0);
      jac.coeffRef(0, 1) =  1.0;
      jac.coeffRef(1, 0) = -1.0;
      jac.coeffRef(1, 1) = 3.0*x(1)*x(1);
    }
  }
};

template <FailType failType>
void run_impl()
{
  using namespace pressio;
  namespace pnonls = pressio::nonlinearsolvers;

  using problem_t  = AdaptedProblem1<failType>;
  using state_t    = typename problem_t::state_type;
  using jacobian_t = typename problem_t::jacobian_type;

  using lin_solver_t = linsol::Solver<linsol::iterative::LSCG, jacobian_t>;
  lin_solver_t linearSolverObj;

  problem_t sys;
  state_t y(2);
  auto nonLinSolver = pnonls::create_newton_solver(sys, linearSolverObj);

  if constexpr (failType == FailType::MaximumIterations) {
    // The max iteration will be hit immediately.
    nonLinSolver.setMaxIterations(1);
  } else if constexpr (failType == FailType::LineSearchStepTooSmall) {
    // First we make sure we are using an Updater that can throw this
    // exception. Then we make the backtrack condition tiny, so alpha
    // keeps getting smaller and smaller until the exception is thrown.
    const auto updateMethod = nonlinearsolvers::Update::BacktrackStrictlyDecreasingObjective;
    nonLinSolver.setUpdateCriterion(updateMethod);
    nonLinSolver.addLineSearchParameter(1e-10);
  } else if constexpr (failType == FailType::LineSearchObjFunctionChangeTooSmall) {
    // Armijo is the only updater that can throw this exception.
    // The rest of the failure should be caused by residualAndJacobian function.
    const auto updateMethod = nonlinearsolvers::Update::Armijo;
    nonLinSolver.setUpdateCriterion(updateMethod);
  } else {
    // Standard values from Problem1 test
    y(0) = 0.001; y(1) = 0.0001;
  }

  try {
    nonLinSolver.solve(y);
  } catch (...) {}

  // Ensure that the file was written out
  std::string fileName = "nonlinear_solver_failed.txt";
  ASSERT_TRUE(std::filesystem::exists(fileName));

  // Read the file to make sure it contains the correct error
  int rank = 0;
#if defined PRESSIO_ENABLE_TPL_MPI
  int flag = 0; MPI_Initialized( &flag );
  if (flag == 1) MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if (rank == 0) {
    std::ifstream file(fileName);
    std::string failMsg; std::getline(file, failMsg);
    ASSERT_NE(failMsg.find(failMsgs.at(failType)), std::string::npos);
    file.close();
    std::filesystem::remove(fileName);
  }
}

TEST(solvers_nonlinear, problem1_residual_failure) {
    run_impl<FailType::ResidualEvaluationFailureUnrecoverable>();
}

TEST(solvers_nonlinear, problem1_residual_has_nans){
  run_impl<FailType::ResidualHasNans>();
}

TEST(solvers_nonlinear, problem1_max_iterations){
  run_impl<FailType::MaximumIterations>();
}

TEST(solvers_nonlinear, problem1_small_step) {
  run_impl<FailType::LineSearchStepTooSmall>();
}

// CANNOT TEST: RootFinder does not support Armijo,
// which is the only updater that would throw this error.
//
// TEST(solvers_nonlinear, problem1_small_obj_ftn_change) {
//   run_impl<FailType::LineSearchObjFunctionChangeTooSmall>();
// }
