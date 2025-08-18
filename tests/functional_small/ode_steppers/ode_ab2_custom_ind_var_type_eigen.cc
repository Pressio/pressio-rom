
#include <gtest/gtest.h>
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "testing_apps.hpp"
#include "custom_independent_variable.hpp"

namespace {

struct AB2MyApp
{
  using independent_variable_type = MyCustomTime;
  using state_type = Eigen::VectorXd;
  using rhs_type = state_type;

  state_type createState() const{ return state_type(3); }

  rhs_type createRhs() const{
    rhs_type R(3);
    return R;
  };

  void rhs(const state_type & /*unused*/,
	   independent_variable_type evaltime,
	   rhs_type & f) const
  {
    f.setConstant(evaltime);
  };
};

struct Collector
{
  void operator()(const ::pressio::ode::StepCount & stepIn,
		  double /*unused*/,
		  const Eigen::VectorXd & y) const
  {
    const auto step = stepIn.get();
    if (step==0){
      EXPECT_DOUBLE_EQ( y(0), 1.);
      EXPECT_DOUBLE_EQ( y(1), 2.);
      EXPECT_DOUBLE_EQ( y(2), 3.);
    }

    if (step==1){
      EXPECT_DOUBLE_EQ( y(0), 1.);
      EXPECT_DOUBLE_EQ( y(1), 2.);
      EXPECT_DOUBLE_EQ( y(2), 3.);
    }

    if (step==2){
      EXPECT_DOUBLE_EQ( y(0), 7.);
      EXPECT_DOUBLE_EQ( y(1), 8.);
      EXPECT_DOUBLE_EQ( y(2), 9.);
    }

    if (step==3){
      EXPECT_DOUBLE_EQ( y(0), 17.);
      EXPECT_DOUBLE_EQ( y(1), 18.);
      EXPECT_DOUBLE_EQ( y(2), 19.);
    }
  }
};
} // end anon namespace

TEST(ode_explicit_steppers, ab2_custom_ind_var)
{
  using namespace pressio;
  using app_t = AB2MyApp;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  y(0) = 1.; y(1) = 2.; y(2) = 3.;

  auto stepperObj = ode::create_ab2_stepper(appObj);

  MyCustomTime t0{0.};
  MyCustomTime dt{2.};
  Collector C;
  auto bp = pressio::ode::steps_fixed_dt(t0, pressio::ode::StepCount(3), dt);
  ode::advance(stepperObj, y, bp, C);
}
