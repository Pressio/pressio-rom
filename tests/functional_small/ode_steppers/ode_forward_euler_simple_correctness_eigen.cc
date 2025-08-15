
#include <gtest/gtest.h>
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "testing_apps.hpp"

#define TEST_ODE_EULER_NUMERICS(y, stepperObj) \
  y(0) = 1.; y(1) = 2.; y(2) = 3.; \
  double dt = 0.1; \
  auto bp = pressio::ode::steps_fixed_dt(0.0, pressio::ode::StepCount(1), dt); \
  ode::advance(stepperObj, y, bp); \
  EXPECT_DOUBLE_EQ( y(0), 1.1); \
  EXPECT_DOUBLE_EQ( y(1), 2.2); \
  EXPECT_DOUBLE_EQ( y(2), 3.3); \
  ode::advance(stepperObj, y, bp); \
  EXPECT_DOUBLE_EQ( y(0), 1.21); \
  EXPECT_DOUBLE_EQ( y(1), 2.42); \
  EXPECT_DOUBLE_EQ( y(2), 3.63); \

TEST(ode_explicit_steppers, euler_system_reference)
{
  using namespace pressio;
  using app_t = ode::testing::AppEigenA;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  auto stepperObj = ode::create_forward_euler_stepper(appObj);
  TEST_ODE_EULER_NUMERICS(y, stepperObj);
}

TEST(ode_explicit_steppers, euler_system_move)
{
  using namespace pressio;
  using app_t = ode::testing::AppEigenA;
  using state_t = typename app_t::state_type;
  app_t appObj;

  state_t y(3);
  auto stepperObj = ode::create_forward_euler_stepper(std::move(appObj));
  TEST_ODE_EULER_NUMERICS(y, stepperObj);
}
