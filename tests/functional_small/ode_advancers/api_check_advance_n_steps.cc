
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

struct AuxClass {
  ScalarType value() { return -8.0; }
  ScalarType value() const { return -9.0; }
};

struct AuxClass2{
  ScalarType value() { return -10.0; }
  ScalarType value()const { return -12.0; }
};

struct Stepper2a{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> tn,
		  pressio::ode::StepCount stepId,
		  pressio::ode::StepSize<independent_variable_type> dt,
		  AuxClass && aux)
  {
    ASSERT_EQ(stepId.get(), 1);
    ASSERT_EQ(tn.get(), 1.0);
    state.fill(2.1 + aux.value() + dt.get());
  }
};

struct Stepper2b{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> tn,
		  pressio::ode::StepCount stepId,
		  pressio::ode::StepSize<independent_variable_type> dt,
		  AuxClass & aux)
  {
    ASSERT_EQ(stepId.get(), 1);
    ASSERT_EQ(tn.get(), 1.0);
    state.fill(2.2 + aux.value() + dt.get());
  }
};

struct Stepper2c{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> tn,
		  pressio::ode::StepCount stepId,
		  pressio::ode::StepSize<independent_variable_type> dt,
		  const AuxClass & aux)
  {
    ASSERT_EQ(stepId.get(), 1);
    ASSERT_EQ(tn.get(), 1.0);
    state.fill(2.3 + aux.value() + dt.get());
  }
};

struct Stepper3a{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> tn,
		  pressio::ode::StepCount stepId,
		  pressio::ode::StepSize<independent_variable_type> dt,
		  AuxClass && aux,
		  AuxClass2 && aux2)
  {
    ASSERT_EQ(stepId.get(), 1);
    ASSERT_EQ(tn.get(), 1.0);
    state.fill(3.1 + aux.value() + aux2.value() + dt.get());
  }
};

struct Stepper3b{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> tn,
		  pressio::ode::StepCount stepId,
		  pressio::ode::StepSize<independent_variable_type> dt,
		  AuxClass && aux,
		  AuxClass2 & aux2)
  {
    ASSERT_EQ(stepId.get(), 1);
    ASSERT_EQ(tn.get(), 1.0);
    state.fill(3.2 + aux.value() + aux2.value() + dt.get());
  }
};

struct Stepper3c{
  using state_type = VectorType;
  using independent_variable_type = ScalarType;
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> tn,
		  pressio::ode::StepCount stepId,
		  pressio::ode::StepSize<independent_variable_type> dt,
		  AuxClass && aux,
		  const AuxClass2 & aux2)
  {
    ASSERT_EQ(stepId.get(), 1);
    ASSERT_EQ(tn.get(), 1.0);
    state.fill(3.3 + aux.value() + aux2.value() + dt.get());
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

#if 0
  {
    Stepper2a stepper;

    advance_n_steps(stepper, odeState, t0, dt, numSteps, AuxClass{});
    check_state_and_reset(odeState, -8.0+2.1+2.0);

    advance_n_steps(stepper, odeState, t0, dtPol, numSteps, AuxClass{});
    check_state_and_reset(odeState, -8.0+2.1+4.0);

    advance_n_steps(stepper, odeState, t0, dt, numSteps, obs, AuxClass{});
    check_state_and_reset(odeState, -8.0+2.1+dt);

    advance_n_steps(stepper, odeState, t0, dtPol, numSteps, obs, AuxClass{});
    check_state_and_reset(odeState, -8.0+2.1+4.0);

    advance_n_steps(stepper, odeState, t0,
		    [](const StepCount & /*unused*/, const StepStartAt<double> & /*unused*/,
		       StepSize<double> & dt){dt = 6.0;},
		    numSteps, AuxClass{});
    check_state_and_reset(odeState, -8.0+2.1+6.0);

    advance_n_steps(stepper, odeState, t0, dt, numSteps,
		    [](const StepCount & /*unused*/, ScalarType /*unused*/,
		       const VectorType & /*unused*/){},
		    AuxClass{});
    check_state_and_reset(odeState, -8.0+2.1+2.0);

    advance_n_steps(stepper, odeState, t0,
		    [](const StepCount & /*unused*/, const StepStartAt<double> & /*unused*/,
		       StepSize<double> & dt){dt=6.0;},
		    numSteps,
		    [](const StepCount & /*unused*/, ScalarType /*unused*/,
		       const VectorType & /*unused*/){},
		    AuxClass{});
    check_state_and_reset(odeState, -8.0+2.1+6.0);
  }

  {
    Stepper2b stepper; AuxClass aux;

    advance_n_steps(stepper, odeState, t0, dt, numSteps, aux);
    check_state_and_reset(odeState, -8.0+2.2+dt);

    advance_n_steps(stepper, odeState, t0, dtPol, numSteps, aux);
    check_state_and_reset(odeState, -8.0+2.2+4.0);

    advance_n_steps(stepper, odeState, t0, dt, numSteps, obs, aux);
    check_state_and_reset(odeState, -8.0+2.2+dt);

    advance_n_steps(stepper, odeState, t0, dtPol, numSteps, obs, aux);
    check_state_and_reset(odeState, -8.0+2.2+4.0);

    advance_n_steps(stepper, odeState, t0,
		    [](const StepCount & /*unused*/, const StepStartAt<double> & /*unused*/,
		       StepSize<double> & dt){dt=5.0;},
		    numSteps, aux);
    check_state_and_reset(odeState, -8.0+2.2+5.0);

    advance_n_steps(stepper, odeState, t0, dt, numSteps,
		    [](const StepCount & /*unused*/, ScalarType /*unused*/, const VectorType & /*unused*/){},
		    aux);
    check_state_and_reset(odeState, -8.0+2.2+dt);

    advance_n_steps(stepper, odeState, t0,
		    [](const StepCount & /*unused*/, const StepStartAt<double> & /*unused*/,
		       StepSize<double> & dt){dt=6.0;},
		    numSteps,
		    [](const StepCount & /*unused*/, ScalarType /*unused*/, const VectorType & /*unused*/){},
		    aux);
    check_state_and_reset(odeState, -8.0+2.2+6.0);
  }

  {
    Stepper2c stepper; const AuxClass aux;

    advance_n_steps(stepper, odeState, t0, dt, numSteps, aux);
    check_state_and_reset(odeState, -9.0+2.3+dt);

    advance_n_steps(stepper, odeState, t0, dtPol, numSteps, aux);
    check_state_and_reset(odeState, -9.0+2.3+4.0);

    advance_n_steps(stepper, odeState, t0, dt, numSteps, obs, aux);
    check_state_and_reset(odeState, -9.0+2.3+dt);

    advance_n_steps(stepper, odeState, t0, dtPol, numSteps, obs, aux);
    check_state_and_reset(odeState, -9.0+2.3+4.0);

    advance_n_steps(stepper, odeState, t0,
		    [](const StepCount & /*unused*/, const StepStartAt<double> & /*unused*/,
		       StepSize<double> & dt){dt=5.0;},
		    numSteps, aux);
    check_state_and_reset(odeState, -9.0+2.3+5.0);

    advance_n_steps(stepper, odeState, t0, dt, numSteps,
		    [](const StepCount & /*unused*/, ScalarType /*unused*/, const VectorType & /*unused*/){},
		    aux);
    check_state_and_reset(odeState, -9.0+2.3+dt);

    advance_n_steps(stepper, odeState, t0,
		    [](const StepCount & /*unused*/, const StepStartAt<double> & /*unused*/,
		       StepSize<double> & dt){dt=6.0;},
		    numSteps,
		    [](const StepCount & /*unused*/, ScalarType /*unused*/, const VectorType & /*unused*/){},
		    aux);
    check_state_and_reset(odeState, -9.0+2.3+6.0);
  }

  {
    Stepper3a stepper; AuxClass2 aux2;

    advance_n_steps(stepper, odeState, t0, dt, numSteps, AuxClass{}, std::move(aux2));
    check_state_and_reset(odeState, -18.0+3.1+2.0);

    advance_n_steps(stepper, odeState, t0, dtPol, numSteps, AuxClass{}, std::move(aux2));
    check_state_and_reset(odeState, -18.0+3.1+4.0);

    advance_n_steps(stepper, odeState, t0, dt, numSteps, obs, AuxClass{}, std::move(aux2));
    check_state_and_reset(odeState, -18.0+3.1+2.0);

    advance_n_steps(stepper, odeState, t0, dtPol, numSteps, obs, AuxClass{}, std::move(aux2));
    check_state_and_reset(odeState, -18.0+3.1+4.0);

    advance_n_steps(stepper, odeState, t0, dt, numSteps,
		    [](const StepCount & /*unused*/, ScalarType /*unused*/, const VectorType & /*unused*/){},
		    AuxClass{}, std::move(aux2));
    check_state_and_reset(odeState, -18.0+3.1+2.0);

    advance_n_steps(stepper, odeState, t0,
		    [](const StepCount & /*unused*/, const StepStartAt<double> & /*unused*/,
		       StepSize<double> & dt){dt=5.0;},
		    numSteps,
		    [](const StepCount & /*unused*/, ScalarType /*unused*/, const VectorType & /*unused*/){},
		    AuxClass{}, std::move(aux2));
    check_state_and_reset(odeState, -18.0+3.1+5.0);
  }

  {
    Stepper3b stepper; AuxClass2 aux2;

    advance_n_steps(stepper, odeState, t0, dt, numSteps, AuxClass{}, aux2);
    check_state_and_reset(odeState, -18.0+3.2+2.0);

    advance_n_steps(stepper, odeState, t0, dtPol, numSteps, AuxClass{}, aux2);
    check_state_and_reset(odeState, -18.0+3.2+4.0);

    advance_n_steps(stepper, odeState, t0, dt, numSteps, obs, AuxClass{}, aux2);
    check_state_and_reset(odeState, -18.0+3.2+2.0);

    advance_n_steps(stepper, odeState, t0, dtPol, numSteps, obs, AuxClass{}, aux2);
    check_state_and_reset(odeState, -18.0+3.2+4.0);

    advance_n_steps(stepper, odeState, t0,
		    [](const StepCount & /*unused*/, const StepStartAt<double> & /*unused*/,
		       StepSize<double> & dt){dt=0.0;},
		    numSteps, AuxClass{}, aux2);
    check_state_and_reset(odeState, -18.0+3.2);
  }

  {
    Stepper3c stepper; const AuxClass2 aux2;

    advance_n_steps(stepper, odeState, t0, dt, numSteps, AuxClass{}, aux2);
    check_state_and_reset(odeState, -20.0+3.3+2.0);

    advance_n_steps(stepper, odeState, t0, dtPol, numSteps, AuxClass{}, aux2);
    check_state_and_reset(odeState, -20.0+3.3+4.0);

    advance_n_steps(stepper, odeState, t0, dt, numSteps, obs, AuxClass{}, aux2);
    check_state_and_reset(odeState, -20.0+3.3+2.0);

    advance_n_steps(stepper, odeState, t0, dtPol, numSteps, obs, AuxClass{}, aux2);
    check_state_and_reset(odeState, -20.0+3.3+4.0);

    advance_n_steps(stepper, odeState, t0,
		    [](const StepCount & /*unused*/, const StepStartAt<double> & /*unused*/,
		       StepSize<double> & dt){dt=0.0;},
		    numSteps, AuxClass{}, aux2);
    check_state_and_reset(odeState, -20.0+3.3);
  }
#endif

}
