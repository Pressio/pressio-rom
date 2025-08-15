
#include <gtest/gtest.h>
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"

struct MyApp1
{
  using independent_variable_type = double;
  using state_type           = Eigen::VectorXd;
  using rhs_type = state_type;
  using mass_matrix_type     = Eigen::MatrixXd;

  state_type createState() const{
    state_type ret(3); ret.setZero();
    return ret;
  }

  rhs_type createRhs() const{
    rhs_type ret(3); ret.setZero();
    return ret;
  };

  mass_matrix_type createMassMatrix() const{
    mass_matrix_type ret(3,3); ret.setZero();
    return ret;
  };

  void massMatrixAndRhs(const state_type & y,
			independent_variable_type evaltime,
			mass_matrix_type & M,
			rhs_type & rhs) const
  {
    for (int i=0; i<y.size(); i++){
      rhs[i] = evaltime + y[i];
    }

    for (int i=0; i<M.rows(); i++){
      for (int j=0; j<M.cols(); j++){
	M(i,j) = evaltime + (double) i;
      }
    }
  };
};

struct FakeLinearSolver1
{
  void solve(const Eigen::MatrixXd & A,
	     Eigen::VectorXd & x,
	     const Eigen::VectorXd & b)
  {
    // we don't need to do anything meaningful here
    // we just need numbers. To make things simple
    // pretend that solving this system boils down
    // to doing x = Ab
    x = A*b;
  }
};

TEST(ode_explicit_steppers, forward_euler_with_mass_matrix)
{
  // this test checkes that things match
  // numbers precomputed manually

  using namespace pressio;
  using app_t   = MyApp1;
  app_t appObj;
  auto stepperObj = ode::create_forward_euler_stepper(appObj);

  FakeLinearSolver1 solver;

  using state_t = typename app_t::state_type;
  {
    state_t y(3);
    y(0) = 1.; y(1) = 2.; y(2) = 3.;
    double dt = 2.;
    auto bp = pressio::ode::steps_fixed_dt(0.0, pressio::ode::StepCount(1), dt);
    ode::advance(stepperObj, y, bp, solver);
    EXPECT_DOUBLE_EQ( y(0), 1.);
    EXPECT_DOUBLE_EQ( y(1), 14.);
    EXPECT_DOUBLE_EQ( y(2), 27.);
  }

  {
    state_t y(3);
    y(0) = 1.; y(1) = 2.; y(2) = 3.;
    double dt = 2.;
    auto bp = pressio::ode::steps_fixed_dt(0.0, pressio::ode::StepCount(2), dt);
    ode::advance(stepperObj, y, bp, solver);
    EXPECT_DOUBLE_EQ( y(0), 193.);
    EXPECT_DOUBLE_EQ( y(1), 302.);
    EXPECT_DOUBLE_EQ( y(2), 411.);
  }

  {
    state_t y(3);
    y(0) = 1.; y(1) = 2.; y(2) = 3.;
    double dt = 2.;
    auto bp = pressio::ode::steps_fixed_dt(0.0, pressio::ode::StepCount(3), dt);    
    ode::advance(stepperObj, y, bp, solver);
    EXPECT_DOUBLE_EQ( y(0), 7344.+ 193.);
    EXPECT_DOUBLE_EQ( y(1), 9180.+ 302.);
    EXPECT_DOUBLE_EQ( y(2), 11016.+411.);
  }
}
