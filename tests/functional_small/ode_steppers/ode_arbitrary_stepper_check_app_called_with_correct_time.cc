
#include <gtest/gtest.h>
#include "pressio/ode_steppers.hpp"
#include "pressio/ode_advancers.hpp"

namespace
{
struct MyApp3
{
  mutable int count_ = 0;
  using independent_variable_type = double;
  using state_type    = Eigen::VectorXd;
  using discrete_residual_type = Eigen::VectorXd;
  using discrete_jacobian_type = Eigen::SparseMatrix<double>;

public:
  state_type createState() const{
    return state_type(3);
  }

  discrete_residual_type createDiscreteResidual() const{
    discrete_residual_type R(3);
    return R;
  }
  discrete_jacobian_type createDiscreteJacobian() const{
    discrete_jacobian_type J(3,3);
    return J;
  }

  template <typename step_t>
  void discreteResidualAndJacobian(const step_t & step,
				   const independent_variable_type & time,
				   const independent_variable_type & dt,
				   discrete_residual_type & /*unused*/,
				   std::optional<discrete_jacobian_type*> /*J*/,
				   const state_type & /*unused*/,
				   const state_type & /*unused*/) const
  {
    ++count_;
    std::cout << count_ << " " << time << " " << dt << std::endl;

    EXPECT_DOUBLE_EQ(dt, 2.5);
    // count_==1: we are doing first step: t_0 = 0 -> t_1
    // so evaluation time should be at t_1
    if(count_==1){
      EXPECT_EQ(step, 1);
      EXPECT_DOUBLE_EQ(time, dt);
    }

    // count_==2: we are doing step: t_1 -> t_2
    // so evaluation time should be at t_2
    if(count_==2){
      EXPECT_EQ(step, 2);
      EXPECT_DOUBLE_EQ(time, 2.*dt);
    }
  }
};

template<typename R_t, typename J_t>
struct MyFakeSolver3
{
  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    R_t R(3);
    J_t J(3,3);
    sys.residualAndJacobian(state, R, &J);
  }
};
}

TEST(ode, mplicit_arbitrary_callWithCorrectTime1)
{
  // fake test to make sure the integrator calls
  // the model with the correct time

  using namespace pressio;
  using app_t = MyApp3;
  using state_t  = typename app_t::state_type;
  app_t appObj;
  state_t y(3); y.setConstant(1.);

  auto stepperObj = ode::create_implicit_stepper<2>(appObj);

  using res_t  = typename app_t::discrete_residual_type;
  using jac_t = typename app_t::discrete_jacobian_type;
  MyFakeSolver3<res_t, jac_t> solver;
  double dt = 2.5;
  auto policy = pressio::ode::steps_fixed_dt(0.0, ode::StepCount(2), dt);
  ode::advance(stepperObj, y, policy, solver);
}

TEST(ode, implicit_arbitrary_callWithCorrectTime2)
{
  using namespace pressio;
  using app_t = MyApp3;
  using state_t  = typename app_t::state_type;
  using res_t  = typename app_t::discrete_residual_type;
  using jac_t = typename app_t::discrete_jacobian_type;
  app_t appObj;
  state_t y(3); y.setConstant(1.);

  auto stepperObj = ode::create_implicit_stepper<2>(appObj);
  MyFakeSolver3<res_t, jac_t> solver;
  double dt = 2.5;
  auto policy = pressio::ode::steps_fixed_dt(0.0, ode::StepCount(2), dt);
  ode::advance(stepperObj, y, policy, solver);
}

TEST(ode, implicit_arbitrary_callWithCorrectTime3)
{
  using namespace pressio;
  using app_t = MyApp3;
  using state_t  = typename app_t::state_type;
  using res_t  = typename app_t::discrete_residual_type;
  using jac_t = typename app_t::discrete_jacobian_type;
  app_t appObj;
  state_t y(3); y.setConstant(1.);

  auto stepperObj = ode::create_implicit_stepper<2>(appObj);
  MyFakeSolver3<res_t, jac_t> solver;
  double dt = 2.5;
  auto policy = pressio::ode::steps_fixed_dt(0.0, ode::StepCount(2), dt);
  ode::advance(stepperObj, y, policy, solver);
}
