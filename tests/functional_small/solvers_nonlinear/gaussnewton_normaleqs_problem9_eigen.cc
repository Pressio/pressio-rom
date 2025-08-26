
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear_gaussnewton.hpp"
#include "./problems/problem9.hpp"

#include <iomanip>

int main()
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  using namespace pressio;
  namespace pnonls = pressio::nonlinearsolvers;

  using problem_t = solvers::test::Problem9<double>;
  using hessian_t = Eigen::MatrixXd;

  problem_t problem;
  auto x = problem.createState();
  x(0) = 1.3;
  x(1) = 6.5e-1;
  x(2) = 6.5e-1;
  x(3) = 7.0e-1;
  x(4) = 6.0e-1;
  x(5) = 3.0;
  x(6) = 5.0;
  x(7) = 7.0;
  x(8) = 2.0;
  x(9) = 4.5;
  x(10) = 5.5;

  using lin_tag      = linearsolvers::direct::HouseholderQR;
  using lin_solver_t = linearsolvers::Solver<lin_tag, hessian_t>;
  lin_solver_t linSolver;

  auto GNSolver = pnonls::create_gauss_newton_solver(problem, linSolver);
  GNSolver.setUpdateCriterion(nonlinearsolvers::Update::BacktrackStrictlyDecreasingObjective);
  GNSolver.setStopTolerance(1e-6);
  GNSolver.solve(problem, x);
  std::cout << std::setprecision(14) << x << " ";
  std::cout << std::endl;

  // this gold solution corresponds to objective = 2.006887e-02
  // which matches the one in http://ftp.mcs.anl.gov/pub/tech_reports/reports/P153.pdf
  // 2.006887e-02 = 0.5 * 0.2003440**2

  Eigen::VectorXd gold(11);
  gold << 1.3099771912462, 0.43155388269712,
    0.63366173059408, 0.59943061295752,
    0.75418346724223, 0.90428784000347,
    1.3658122743845,  4.8236981348807,
    2.3986848313274,  4.5688746007188,
    5.6753414536432;

  if (gold.isApprox(x)){
    std::cout << "PASSED" << std::endl;
  } else{
    std::cout << "FAILED" << std::endl;
  }

  PRESSIOLOG_FINALIZE();
  return 0;
}
