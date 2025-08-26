
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear_gaussnewton.hpp"
#include "./problems/problem3.hpp"

#include <iomanip>

template<class scalar_t>
struct IdentityWeigher{

  template<class T>
  void operator()(const Eigen::MatrixBase<T> & operand,
		  Eigen::MatrixBase<T> & result) const
  {
    // since this is an identity, just do a deep copy
    result = operand;
  }
};

template <typename problem_t, typename state_t, typename solver>
void testC1(std::string & sentinel,
            problem_t & problem,
            state_t & x,
            solver & GNSolver)
{
  GNSolver.setStopTolerance(1e-8);
  GNSolver.solve(problem, x);
  if ( (std::abs(x(0) - 2.4173449278229) > 1e-7 )or
       (std::abs(x(1) - 0.26464986197941) > 1e-7) ){
    sentinel = "FAILED";
  }
}

template <typename problem_t, typename state_t, typename solver>
void testC2(std::string & sentinel,
            problem_t & problem,
            state_t & x,
            solver & GNSolver)
{
  using namespace pressio::nonlinearsolvers;
  auto criterion = Stop::WhenAbsolutel2NormOfGradientBelowTolerance;
  // 1e3 is chosen to test the convergence condition
  GNSolver.setStopCriterion(criterion);
  GNSolver.setStopTolerance(1e3);
  GNSolver.solve(problem, x);

  if ( (std::abs(x(0) - 2.4153884777105201) > 1e-11 )or
       (std::abs(x(1) - 0.26879930127342189) > 1e-11) ){
    sentinel = "FAILED";
  }
}

template <typename problem_t, typename state_t, typename solver>
void testC3(std::string & sentinel,
            problem_t & problem,
            state_t & x,
            solver & GNSolver)
{
  using namespace pressio::nonlinearsolvers;
  auto criterion = Stop::WhenRelativel2NormOfGradientBelowTolerance;
  GNSolver.setStopCriterion(criterion);
  GNSolver.setStopTolerance(1e-5);
  GNSolver.solve(problem, x);
  if ( (std::abs(x(0) - 2.41728158794844) > 1e-11 )or
       (std::abs(x(1) - 0.26465375115797) > 1e-11) ){
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
  using mat_t   = typename problem_t::jacobian_type;
  using hessian_t = mat_t;

  problem_t problem;
  state_t x(2); x(0) = 2.0; x(1) = 0.25;

  // linear solver type
  using solver_tag	= linsol::iterative::LSCG;
  using linear_solver_t = linsol::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolver;

  auto GNSolver = pnonls::create_gauss_newton_solver(problem, linSolver, IdentityWeigher<double>{});

  x(0) = 2.0; x(1) = 0.25;
  testC1(sentinel, problem, x, GNSolver);
  std::cout << "\n" << std::endl;

  x(0) = 2.0; x(1) = 0.25;
  testC2(sentinel, problem, x, GNSolver);
  std::cout << "\n" << std::endl;

  x(0) = 2.0; x(1) = 0.25;
  testC3(sentinel, problem, x, GNSolver);
  std::cout << "\n" << std::endl;

  std::cout << sentinel << std::endl;

  std::cout << std::setprecision(14) << x << std::endl;
  PRESSIOLOG_FINALIZE();
  return 0;
}
