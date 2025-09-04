
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear_gaussnewton.hpp"
#include "./problems/problem5.hpp"

#include <iomanip>

int main()
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::info);

  using namespace pressio;

  Eigen::Vector4d state;

  using problem_t = solvers::test::Problem5a<double>;
  using state_t   = typename problem_t::state_type;
  using hessian_t = Eigen::MatrixXd;

  problem_t problem;
  state_t x(4);
  x(0) = -0.05; x(1) = 1.1; x(2) = 1.2; x(3) = 1.5;

  using lin_tag      = linsol::direct::HouseholderQR;
  using lin_solver_t = linsol::Solver<lin_tag, hessian_t>;
  lin_solver_t linSolver;

  auto GNSolver = nlsol::create_gauss_newton_solver(problem, linSolver);
  GNSolver.setStopTolerance(1e-5);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << x << std::endl;

  std::vector<double> gold = {1.00000001567414e+00,
			      9.99999999124769e-01,
			      9.99999996519930e-01,
			      9.99999988898883e-01};

  std::string sentinel = "PASSED";
  const auto e1 = std::abs(x(0) - gold[0]);
  const auto e2 = std::abs(x(1) - gold[1]);
  const auto e3 = std::abs(x(2) - gold[2]);
  const auto e4 = std::abs(x(3) - gold[3]);
  if (e1>1e-6 or e2>1e-6  or e3>1e-6 or e4>1e-6){
    sentinel = "FAILED";
  }
  std::cout << sentinel << std::endl;

  PRESSIOLOG_FINALIZE();
  return 0;
}
