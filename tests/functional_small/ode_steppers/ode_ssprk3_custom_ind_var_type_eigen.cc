
#include <gtest/gtest.h>
#include "pressio/ode_steppers.hpp"
#include "pressio/ode_advancers.hpp"
#include "custom_independent_variable.hpp"

struct AppForSSPRK3
{
  using independent_variable_type = MyCustomTime;
  using state_type = Eigen::VectorXd;
  using rhs_type = state_type;

  state_type createState() const{
    return rhs_type(3);
  };

  rhs_type createRhs() const{
    return rhs_type(3);
  };

  void rhs(const state_type & y,
	   const independent_variable_type evaltime,
	   rhs_type & rhs) const
  {
    auto sz = y.size();
    for (decltype(sz) i=0; i<sz; i++){
      rhs[i] = y[i] + evaltime;
    }
  };
};

TEST(ode_explicit_steppers, ssprk3_custom_ind_var)
{
  using namespace pressio;
  using app_t = AppForSSPRK3;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  y(0) = 1.; y(1) = 2.; y(2) = 3.;
  auto stepperObj = ode::create_ssprk3_stepper(appObj);

  MyCustomTime t0{0.};
  MyCustomTime dt{2.};
  auto bp = pressio::ode::steps_fixed_dt(t0, pressio::ode::StepCount(1), dt);
  ode::advance(stepperObj, y, bp);
  EXPECT_DOUBLE_EQ( y(0), 29./3.);
  EXPECT_DOUBLE_EQ( y(1), 48./3.);
  EXPECT_DOUBLE_EQ( y(2), 67./3.);
}
