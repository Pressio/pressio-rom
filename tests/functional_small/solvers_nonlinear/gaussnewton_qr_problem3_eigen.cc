
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear_gaussnewton.hpp"
#include "./problems/problem3.hpp"

#include <iomanip>

template <typename state_t, typename solver>
void testC1(std::string & sentinel,
            state_t & x,
            solver & GNSolver)
{
  GNSolver.setStopTolerance(1e-8);
  GNSolver.solve(x);
  std::cout << std::setprecision(16) << x << std::endl;
  if ( (std::abs(x(0) - 2.4173449278229) > 1e-7 )or
       (std::abs(x(1) - 0.26464986197941) > 1e-7) ){
    sentinel = "FAILED";
  }
}

template <typename state_t, typename solver>
void testC2(std::string & sentinel,
            state_t & x,
            solver & GNSolver)
{
  GNSolver.setMaxIterations(2);
  GNSolver.solve(x);
  std::cout << std::setprecision(16) << x << std::endl;
  if ( (std::abs(x(0) - 2.415361667771343) > 1e-9 )or
       (std::abs(x(1) - 0.2648293802571118) > 1e-9) ){
    sentinel = "FAILED";
  }
}

int main()
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::sparse);

  std::string sentinel = "PASSED";
  using namespace pressio;
  namespace pnonls = pressio::nonlinearsolvers;
  
  using problem_t   = solvers::test::Problem3<double>;
  using state_t	    = typename problem_t::state_type;

  problem_t problem;
  state_t x(2); x(0) = 2.0; x(1) = 0.25;

  using mat_type  = typename problem_t::jacobian_type;
  using qr_solver_type = qr::QRSolver<mat_type, qr::Householder>;
  qr_solver_type qrSolver;

  auto GNSolver = pnonls::experimental::create_gauss_newton_qr_solver(problem, qrSolver);

  x(0) = 2.0; x(1) = 0.25;
  testC1(sentinel, x, GNSolver);
  std::cout << "\n" << std::endl;

  x(0) = 2.0; x(1) = 0.25;
  testC2(sentinel, x, GNSolver);
  std::cout << "\n" << std::endl;

  std::cout << sentinel << std::endl;

  std::cout << std::setprecision(14) << x << std::endl;
  PRESSIOLOG_FINALIZE();
  return 0;
}
