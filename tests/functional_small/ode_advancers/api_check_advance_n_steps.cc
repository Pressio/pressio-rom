
#include <gtest/gtest.h>
#include "pressio/ode_advancers.hpp"

using ScalarType = double;
using VectorType = std::array<ScalarType,3>;

struct Stepper1{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> tn,
		  pressio::ode::StepCount stepId,
		  pressio::ode::StepSize<independent_variable_type> dt)
  {
    ASSERT_EQ(stepId.get(), 1);
    ASSERT_EQ(tn.get(), 1.0);
    state.fill(1. + dt.get());
  }
};

struct Observer{
  void operator()(const pressio::ode::StepCount & /*unused*/,
		  const ScalarType & /*unused*/,
		  const VectorType & /*unused*/) const{}
};

struct DtSetter{
  void operator()(const pressio::ode::StepCount & /*unused*/,
                  const pressio::ode::StepStartAt<double> & /*unused*/,
                  pressio::ode::StepSize<double> & dt) const
  {
    dt = 4.0;
  }
};

void check_state_and_reset(VectorType & a, ScalarType value){
  EXPECT_DOUBLE_EQ(a[0], value);
  EXPECT_DOUBLE_EQ(a[1], value);
  EXPECT_DOUBLE_EQ(a[2], value);
  a.fill(0.);
}

TEST(ode, test)
{
  VectorType odeState; odeState.fill(0.);
  ScalarType dt = 2.0;
  ScalarType t0 = 1.0;
  const auto numSteps = pressio::ode::StepCount(1);
  Observer obs;
  DtSetter dtPol;

  using namespace pressio::ode;

  {
    Stepper1 stepper;

    auto bp = pressio::ode::steps_fixed_dt(t0, numSteps, dt);
    advance(stepper, odeState, bp);
    check_state_and_reset(odeState, 3.0);

    auto bp2 = pressio::ode::steps(t0, numSteps, dtPol);
    advance(stepper, odeState, bp2);
    check_state_and_reset(odeState, 5.0);

    auto bp3 = pressio::ode::steps_fixed_dt(t0, numSteps, dt);
    advance(stepper, odeState, bp3, obs);
    check_state_and_reset(odeState, 3.0);

    auto bp4 = pressio::ode::steps(t0, numSteps, dtPol);
    advance(stepper, odeState, bp4, obs);
    check_state_and_reset(odeState, 5.0);

    auto bp5 = pressio::ode::steps(t0, numSteps,
				   [](const StepCount & /*unused*/,
				      const StepStartAt<double> & /*unused*/,
				      StepSize<double> & dt){
				     dt = 5.0; });
    advance(stepper, odeState, bp5);
    check_state_and_reset(odeState, 6.0);


    auto bp6 = pressio::ode::steps(t0, numSteps,
				   [](const StepCount & /*unused*/,
				      const StepStartAt<double> & /*unused*/,
				      StepSize<double> & dt){
				     dt = 5.0; });
    advance(stepper, odeState, bp6,
	    [](const StepCount & /*unused*/, ScalarType /*unused*/, const VectorType & /*unused*/){});
    check_state_and_reset(odeState, 6.0);
  }
}
