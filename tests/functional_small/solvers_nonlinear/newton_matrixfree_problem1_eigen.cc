
#include <gtest/gtest.h>
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear_newton.hpp"

struct Problem1MatrixFree
{
  using state_type = Eigen::VectorXd;
  using residual_type = state_type;

private:
  mutable Eigen::SparseMatrix<double> J_;

public:
  state_type createState() const {
    state_type a(2);
    a.setZero();
    return a;
  }

  residual_type createResidual() const {
    residual_type a(2);
    a.setZero();
    return a;
  }

  template<class OperandT, class ResultT>
  void applyJacobian(state_type const & x,
		     OperandT const & operand,
		     ResultT & out) const
  {
    // here we compute J but this is just for convenience here,
    // we are only interested in the jacobian action
    J_.resize(2,2);
    J_.coeffRef(0, 0) = 3.0*x(0)*x(0);
    J_.coeffRef(0, 1) =  1.0;
    J_.coeffRef(1, 0) = -1.0;
    J_.coeffRef(1, 1) = 3.0*x(1)*x(1);

    out = J_ * operand;
  }

  void residual(const state_type& x,
		residual_type& res) const
  {
    res(0) =  x(0)*x(0)*x(0) + x(1) - 1.0;
    res(1) = -x(0) + x(1)*x(1)*x(1) + 1.0;
  }
};

TEST(solvers_nonlinear, problem1MatrixFree)
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  using namespace pressio;
  namespace pnonls = pressio::nonlinearsolvers;

  using problem_t  = Problem1MatrixFree;
  using state_t    = typename problem_t::state_type;

  problem_t P;
  state_t y(2);
  using tag = pressio::linearsolvers::iterative::GMRES;
  auto nonLinSolver = pnonls::experimental::create_matrixfree_newtonkrylov_solver<tag>(P);

  y(0) = 0.001;
  y(1) = 0.0001;

  nonLinSolver.solve(y);
  const auto e1 = std::abs(y(0) - (1.));
  const auto e2 = std::abs(y(1) - (0.));
  std::cout << y << std::endl;
  ASSERT_TRUE(e1<1e-8);
  ASSERT_TRUE(e2<1e-8);

  PRESSIOLOG_FINALIZE();
}
