
#include <gtest/gtest.h>
#include "pressio/ode_steppers.hpp"
#include "pressio/ode_advancers.hpp"

namespace
{

struct MyApp1
{
  using independent_variable_type   = double;
  using state_type    = Eigen::VectorXd;
  using rhs_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

public:
  state_type createState() const{
    return state_type(3);
  }

  rhs_type createRhs() const{
    rhs_type f(3);
    return f;
  }

  jacobian_type createJacobian() const{
    jacobian_type J(3,3);
    // ensure that the diagonal elements exist
    for (int i=0; i<J.innerSize(); i++) {
      J.coeffRef(i, i) = 0;
    }
    return J;
  }

  void rhsAndJacobian(const state_type & /*unused*/,
		      const independent_variable_type& evaltime,
		      rhs_type & f,
		      std::optional<jacobian_type*> /*J*/) const
  {
    std::cout << "f: t=" << evaltime << "\n";
    f[0] = evaltime+1.;
    f[1] = evaltime+2.;
    f[2] = evaltime+3.;
  }
};

struct MyFakeSolver1
{
  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    auto R = sys.createResidual();
    auto J = sys.createJacobian();
    sys.residualAndJacobian(state, R, std::optional<decltype(J)*>(&J));

    // here we should have:
    // R = y_n - y_n-1 - dt*f()
    // computedd at t_1 = dt = 1.1
    // where y_n = y_n-1 = [1.5 1.5 1.5]
    // f = [2.1 3.1 4.1]
    //
    // so R = -1.1*[2.1 3.1 4.1] = -[2.31 3.41 4.51]
    std::cout << R << std::endl;
    EXPECT_DOUBLE_EQ(R(0), -2.31);
    EXPECT_DOUBLE_EQ(R(1), -3.41);
    EXPECT_DOUBLE_EQ(R(2), -4.51);
  }
};

struct MyApp2
{
  using independent_variable_type   = double;
  using state_type    = Eigen::VectorXd;
  using rhs_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

public:
  state_type createState() const{
    return state_type(3);
  }

  rhs_type createRhs() const{
    rhs_type f(3);
    return f;
  }

  jacobian_type createJacobian() const{
    jacobian_type J(3,3);
    // ensure that the diagonal elements exist
    for (int i=0; i<J.innerSize(); i++) {
      J.coeffRef(i, i) = 0;
    }
    return J;
  }

  void rhsAndJacobian(const state_type & /*unused*/,
		      const independent_variable_type& evaltime,
		      rhs_type & f,
		      std::optional<jacobian_type*> /*J*/) const
  {
    std::cout << "f: t=" << evaltime << "\n";
    f[0] = evaltime+1.;
    f[1] = evaltime+2.;
    f[2] = evaltime+3.;
  }
};

struct MyFakeSolver2
{
  int count_ = {};

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    ++count_;

    auto R = sys.createResidual();
    auto J = sys.createJacobian();
    sys.residualAndJacobian(state, R, std::optional<decltype(J)*>(&J));

    if (count_==1){
      // this is called from auxiliary stepper which is bdf1

      // here we should have:
      // R = y_n - y_n-1 - dt*f()
      // computedd at t_1 = dt = 1.1
      // where y_n = y_n-1 = [1.5 1.5 1.5]
      // f = [2.1 3.1 4.1]
      EXPECT_DOUBLE_EQ(R(0), 1.5-1.5-1.1*2.1);
      EXPECT_DOUBLE_EQ(R(1), 1.5-1.5-1.1*3.1);
      EXPECT_DOUBLE_EQ(R(2), 1.5-1.5-1.1*4.1);
    }
    if (count_==2)
    {
      // this is called from bdf2
      // here we should have:
      // R = y_n+1 - (4/3)y_n + (1/3)*y_n-1 - (2/3)*dt*f()
      // computedd at t_2 = dt = 2.2
      // where y_n+1 = y_n = y_n-1 = [1.5 1.5 1.5]
      // f = [3.2 4.2 5.2]
      //
      std::cout << state << std::endl;
      std::cout << R << std::endl;
      EXPECT_DOUBLE_EQ(R(0), 1.5-(4./3.)*1.5+(1./3.)*1.5 - (2./3.)*1.1*3.2);
      EXPECT_DOUBLE_EQ(R(1), 1.5-(4./3.)*1.5+(1./3.)*1.5 - (2./3.)*1.1*4.2);
      EXPECT_DOUBLE_EQ(R(2), 1.5-(4./3.)*1.5+(1./3.)*1.5 - (2./3.)*1.1*5.2);
    }
  }
};

} // end anonym namespace

TEST(ode, implicit_euler_appVelocityCalledWithCorrectTime)
{
  // fake test to make sure the app velocity
  // for bdf1 is called with the time at the next time step
  // so here we set the velocity to be time+something
  // so that when the "fake" solver is called we can verify
  // that the residual has the correct numbers

  using namespace pressio;
  using app_t = MyApp1;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3); y.setConstant(1.5);

  auto stepperObj = ode::create_bdf1_stepper(appObj);
  MyFakeSolver1 solver;
  double dt = 1.1;
  auto policy = ode::steps_fixed_dt(0.0, ode::StepCount(1), dt);
  ode::advance(stepperObj, y, policy, solver);
}

TEST(ode, implicit_bdf2_appVelocityCalledWithCorrectTime)
{
  // fake test to make sure the app velocity
  // for bdf1 is called with the time at the next time step
  // so here we set the velocity to be time+something
  // so that when the "fake" solver is called we can verify
  // that the residual has the correct numbers

  using namespace pressio;
  using app_t = MyApp2;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3); y.setConstant(1.5);

  auto stepperObj = ode::create_bdf2_stepper(appObj);
  MyFakeSolver2 solver;
  double dt = 1.1;
  auto policy = ode::steps_fixed_dt(0.0, ode::StepCount(2), dt);  
  ode::advance(stepperObj, y, policy, solver);
}
