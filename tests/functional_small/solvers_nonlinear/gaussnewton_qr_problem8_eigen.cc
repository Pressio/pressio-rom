
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear_gaussnewton.hpp"
#include "./problems/problem8.hpp"

#include <iomanip>

template <typename state_t, typename solver>
void testC1(std::string & sentinel,
            state_t & x,
            solver & GNSolver)
{
  GNSolver.setStopTolerance(1e-8);
  GNSolver.solve(x);
  if ( (std::abs(x(0) - 5.6096364710E-03) > 2.e-4) or
       (std::abs(x(1) - 6.1813463463E+03) > 3.e+1) or
       (std::abs(x(2) - 3.4522363462E+02) > 8.e-1) ){
    sentinel = "FAILED";
  }
}

template <typename state_t, typename solver>
void testC2(std::string & sentinel,
            state_t & x,
            solver & GNSolver)
{
  // setting 2 max iters so the solver will exit after 20 steps
  GNSolver.setMaxIterations(20);
  GNSolver.solve(x);
  std::cout << std::setprecision(16) << *x.data() << std::endl;
  if ( (std::abs(x(0) - 5.6096364710E-03) > 2.e-4) or
       (std::abs(x(1) - 6.1813463463E+03) > 3.e+1) or
       (std::abs(x(2) - 3.4522363462E+02) > 8.e-1) ){
    sentinel = "FAILED";
  }
}

int main()
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::sparse);

  std::string sentinel= "PASSED";
  using namespace pressio;

  using problem_t = solvers::test::Problem8<double>;
  using state_t = typename problem_t::state_type;
  problem_t problem;
  state_t x(3);
  x(0) = 0.02; x(1) = 4000.0; x(2) = 250.0;

  using mat_type  = typename problem_t::jacobian_type;
  using qr_solver_type = qr::QRSolver<mat_type, qr::Householder>;
  qr_solver_type qrSolver;

  auto GNSolver = nlsol::experimental::create_gauss_newton_qr_solver(problem, qrSolver);

  x(0) = 0.02; x(1) = 4000.0; x(2) = 250.0;
  testC1(sentinel, x, GNSolver);
  std::cout << "\n" << std::endl;

  x(0) = 0.02; x(1) = 4000.0; x(2) = 250.0;
  testC2(sentinel, x, GNSolver);
  std::cout << "\n" << std::endl;

  std::cout << sentinel << std::endl;
  PRESSIOLOG_FINALIZE();
  return 0;
}
