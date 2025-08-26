
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear_gaussnewton.hpp"
#include "./problems/problem3.hpp"

#include <iomanip>

template<class scalar_t>
struct RangeWeigher{

  Eigen::Matrix<scalar_t, -1, -1> operator_w;

  RangeWeigher() {
    operator_w.resize(leadingDim(), 8);
    for (int i = 0; i < leadingDim(); ++i) {
      for (int j = 0; j < 8; ++ j) {
        operator_w(i, j) = i;
      }
    }
  }

  int leadingDim() const { return 4; }

  void operator()(const Eigen::Matrix<scalar_t, -1, 1> & operand,
		  Eigen::Matrix<scalar_t, -1, 1> & result) const
  {
    pressio::ops::product(
        pressio::nontranspose(),
        1., operator_w, operand,
        0., result
    );
  }

  void operator()(const Eigen::Matrix<scalar_t, -1, -1> & operand,
		  Eigen::Matrix<scalar_t, -1, -1> & result) const
  {
    pressio::ops::product(
        ::pressio::nontranspose(), ::pressio::nontranspose(),
        1., operator_w, operand,
        0., result
    );
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
  if ( (std::abs(x(0) - 2.2812490140064) > 1e-7 )or
       (std::abs(x(1) - 0.27503679698785) > 1e-7) ){
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

  if ( (std::abs(x(0) - 2.2817056308917141) > 1e-11 )or
       (std::abs(x(1) - 0.27507107642086648) > 1e-11) ){
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
  if ( (std::abs(x(0) - 2.2812490853208414) > 1e-11 )or
       (std::abs(x(1) - 0.27503680234281269) > 1e-11) ){
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
  using tag_t = pnonls::impl::CompactWeightedGaussNewtonNormalEqTag;
  using solver_tag	= linearsolvers::iterative::LSCG;
  using linear_solver_t = linearsolvers::Solver<solver_tag, hessian_t>;
  linear_solver_t linSolver;
  RangeWeigher<double> W{};
  auto GNSolver = pnonls::create_gauss_newton_solver(problem, linSolver, W, tag_t{});

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
