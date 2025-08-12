
#include <gtest/gtest.h>
#include "pressio/ode_advancers.hpp"

using ScalarType = double;
using VectorType = std::vector<ScalarType>;

struct Stepper1
{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;

  void operator()(state_type & odeState,
		  pressio::ode::StepStartAt<independent_variable_type> /*unused*/,
		  pressio::ode::StepCount /*unused*/,
		  pressio::ode::StepSize<independent_variable_type> dt)
  {
    for (std::size_t i=0; i<odeState.size(); i++){
      odeState[i] += dt.get();
    }
  }
};

struct DtSetter1{
  void operator()(pressio::ode::StepCount currStep,
                  pressio::ode::StepStartAt<double> /*unused*/,
                  pressio::ode::StepSize<double> & dt) const
  {
   if(currStep.get()==1) dt = 2.;
   if(currStep.get()==2) dt = 0.5;
   if(currStep.get()==3) dt = 1.;
   if(currStep.get()==4) dt = 1.2;
  }
};

TEST(ode, advance_n_steps_varying_dt_stepper)
{
  VectorType odeState(5);
  std::for_each(odeState.begin(), odeState.end(),
		[](ScalarType & val){val = 0.0; });

  Stepper1 stepper;
  ScalarType t0 = 0.;
  DtSetter1 dtSetter;

  auto bp = pressio::ode::steps(t0, pressio::ode::StepCount(4), dtSetter);
  pressio::ode::advance(stepper, odeState, bp);
  std::for_each(odeState.begin(), odeState.end(),
		[](const ScalarType & val){ EXPECT_DOUBLE_EQ(val, 4.7);});
}
