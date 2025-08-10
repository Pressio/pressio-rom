
#include <gtest/gtest.h>
#include "pressio/ode_concepts.hpp"

namespace{
  using namespace pressio::ode;

  struct Stepper1{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    void operator()(state_type &	    /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepCount		    /*unused*/,
		    StepSize<independent_variable_type> /*unused*/){}
  };

  struct Stepper2{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    // wrong order
    void operator()(state_type &	    /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepSize<independent_variable_type> /*unused*/,
		    StepCount		    /*unused*/){}
  };

  struct Stepper3{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    // missing ref here
    void operator()(state_type              /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepCount               /*unused*/,
		    StepSize<independent_variable_type> /*unused*/){}
  };

  struct Stepper4{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    // rvalue ref here
    void operator()(state_type &&           /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepCount               /*unused*/,
		    StepSize<independent_variable_type> /*unused*/){}
  };


  struct AuxThing1{};
  struct AuxThing2{};

  struct VarStepper1{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    void operator()(state_type &            /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepCount               /*unused*/,
		    StepSize<independent_variable_type> /*unused*/,
		    AuxThing1               /*unused*/){}
  };

  struct VarStepper2{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    // wrong order
    void operator()(state_type &            /*unused*/,
		    StepCount               /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepSize<independent_variable_type> /*unused*/,
		    AuxThing1               /*unused*/){}
  };

  struct VarStepper3{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    // missing ref for state
    void operator()(state_type              /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepCount               /*unused*/,
		    StepSize<independent_variable_type> /*unused*/,
		    AuxThing1               /*unused*/){}
  };

  struct VarStepper4{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    // rvalue ref for state
    void operator()(state_type &&           /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepCount               /*unused*/,
		    StepSize<independent_variable_type> /*unused*/,
		    AuxThing1               /*unused*/){}
  };


  struct VarStepper5{
    using independent_variable_type = double;
    using state_type = std::vector<float>;

    void operator()(state_type &            /*unused*/,
		    StepStartAt<independent_variable_type>    /*unused*/,
		    StepCount               /*unused*/,
		    StepSize<independent_variable_type> /*unused*/,
		    AuxThing1               /*unused*/,
		    AuxThing2               /*unused*/)
    {}
  };

} //end anonym namespace

TEST(ode, concepts__stepper)
{
  using namespace pressio::ode;

  static_assert( ExplicitStepper<Stepper1>::value, "");
  static_assert(!ExplicitStepper<Stepper2>::value, "");
}

TEST(ode, concepts_variadic_stepper)
{
  using namespace pressio::ode;

  static_assert(!ExplicitStepper<VarStepper1>::value, "");
  static_assert(!ExplicitStepper<VarStepper2>::value, "");
  static_assert(!ExplicitStepper<VarStepper3>::value, "");
  static_assert(!ExplicitStepper<VarStepper4>::value, "");
  static_assert(!ExplicitStepper<VarStepper5>::value, "");

  // static_assert(ImplicitStepper<void, VarStepper1, AuxThing1>::value, "");
  // static_assert(!ImplicitStepper<void, VarStepper2, AuxThing1>::value, "");
  // static_assert(!ImplicitStepper<void, VarStepper5, AuxThing1>::value, "");
  // static_assert(ImplicitStepper<void, VarStepper5, AuxThing1, AuxThing2>::value, "");
}
