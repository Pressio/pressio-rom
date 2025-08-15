
#include <gtest/gtest.h>
#include "pressio/ode_advancers.hpp"

using ScalarType = double;
using VectorType = std::vector<ScalarType>;

struct Stepper
{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;

  void operator()(state_type & /*unused*/,
		  pressio::ode::StepStartAt<independent_variable_type> currTime,
		  pressio::ode::StepCount step,
		  pressio::ode::StepSize<independent_variable_type> dt)
  {
    EXPECT_DOUBLE_EQ(dt.get(), 2.);
    if (step.get()==1){
      EXPECT_DOUBLE_EQ(currTime.get(), 0.);
    }
    if (step.get()==2){
      EXPECT_DOUBLE_EQ(currTime.get(), 2.);
    }
    if (step.get()==3){
      EXPECT_DOUBLE_EQ(currTime.get(), 4.);
    }
    if (step.get()==4){
      EXPECT_DOUBLE_EQ(currTime.get(), 6.);
    }
    if (step.get()==5){
      EXPECT_DOUBLE_EQ(currTime.get(), 8.);
    }
  }
};

TEST(ode, advance_n_steps_check_step_and_time)
{
  VectorType odeState(5);
  Stepper stepper;
  ScalarType dt = 2.;
  ScalarType t0 = 0.;
  auto bp = pressio::ode::steps_fixed_dt(t0, pressio::ode::StepCount(5), dt);
  pressio::ode::advance(stepper, odeState, bp);
}
