
#include <gtest/gtest.h>
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"

struct MyApp
{
  using independent_variable_type = double;
  using state_type           = Eigen::VectorXd;
  using rhs_type = state_type;

  mutable int count1 = 0;
  const std::map<int, Eigen::VectorXd> & rhs_;

  MyApp(const std::map<int, Eigen::VectorXd> & rhs)
    : rhs_(rhs){}

  state_type createState() const{
    state_type ret(3); ret.setZero();
    return ret;
  }

  rhs_type createRhs() const{
    rhs_type ret(3); ret.setZero();
    return ret;
  };

  void rhs(const state_type & y,
	   independent_variable_type evaltime,
	   rhs_type & rhs) const
  {
    rhs = rhs_.at(count1++);
  };
};

struct StateObserver
{
  template<class IndepVarType>
  void operator()(pressio::ode::StepCount /**/,
		  IndepVarType /**/,
		  const Eigen::VectorXd & /**/) const
  {}
};

template<int method_switch>
struct RhsObserver
{
  const double dt_ = {};
  const std::map<int, Eigen::VectorXd> & rhs_;
  RhsObserver(double dt, const std::map<int, Eigen::VectorXd> & rhs)
    : dt_(dt), rhs_(rhs){}

  template<class IndepVarType, int _method_switch = method_switch>
  typename std::enable_if<_method_switch == 0>::type
  operator()(pressio::ode::StepCount stepIn,
	     pressio::ode::IntermediateStepCount imStepIn,
	     IndepVarType timeIn,
	     const Eigen::VectorXd & rhsIn) const
  {
    //
    // euler
    //

    // does not have intermediate steps
    EXPECT_TRUE(imStepIn.get() == 0);

    // the -1 is because step starts from 1
    EXPECT_TRUE( rhsIn.isApprox(rhs_.at(stepIn.get()-1)) );

    const IndepVarType expectedTime = (stepIn.get()-1)*dt_;
    EXPECT_TRUE( std::abs(expectedTime - timeIn) < 1e-13 );
  }

  template<class IndepVarType, int _method_switch = method_switch>
  typename std::enable_if<_method_switch == 1>::type
  operator()(pressio::ode::StepCount stepIn,
	     pressio::ode::IntermediateStepCount imStepIn,
	     IndepVarType timeIn,
	     const Eigen::VectorXd & rhsIn) const
  {
    //
    // rk4
    //

    // has intermediate steps
    EXPECT_TRUE(imStepIn.get() == 0 or
		imStepIn.get() == 1 or
		imStepIn.get() == 2 or
		imStepIn.get() == 3);

    // the -1 is because step starts from 1
    const int index = (stepIn.get()-1) * 4 + imStepIn.get();
    EXPECT_TRUE( rhsIn.isApprox(rhs_.at(index)) );

    if (imStepIn.get() == 0){
      const IndepVarType expectedTime = (stepIn.get()-1)*dt_;
      EXPECT_TRUE( std::abs(expectedTime - timeIn) < 1e-13 );
    }
    if (imStepIn.get() == 1){
      const IndepVarType expectedTime = (stepIn.get()-1)*dt_ + 0.5*dt_;
      EXPECT_TRUE( std::abs(expectedTime - timeIn) < 1e-13 );
    }

    if (imStepIn.get() == 2){
      const IndepVarType expectedTime = (stepIn.get()-1)*dt_ + 0.5*dt_;
      EXPECT_TRUE( std::abs(expectedTime - timeIn) < 1e-13 );
    }
    if (imStepIn.get() == 3){
      const IndepVarType expectedTime = (stepIn.get()-1)*dt_ + dt_;
      EXPECT_TRUE( std::abs(expectedTime - timeIn) < 1e-13 );
    }
  }

  template<class IndepVarType, int _method_switch = method_switch>
  typename std::enable_if<_method_switch == 2>::type
  operator()(pressio::ode::StepCount stepIn,
	     pressio::ode::IntermediateStepCount imStepIn,
	     IndepVarType timeIn,
	     const Eigen::VectorXd & rhsIn) const
  {
    //
    // ab2
    //

    // has NO intermediate steps
    EXPECT_TRUE(imStepIn.get() == 0);

    // the -1 is because step starts from 1
    EXPECT_TRUE( rhsIn.isApprox(rhs_.at(stepIn.get()-1)) );

    const IndepVarType expectedTime = (stepIn.get()-1)*dt_;
    EXPECT_TRUE( std::abs(expectedTime - timeIn) < 1e-13 );
  }

  template<class IndepVarType, int _method_switch = method_switch>
  typename std::enable_if<_method_switch == 3>::type
  operator()(pressio::ode::StepCount stepIn,
	     pressio::ode::IntermediateStepCount imStepIn,
	     IndepVarType timeIn,
	     const Eigen::VectorXd & rhsIn) const
  {
    //
    // rk3
    //

    // has intermediate steps
    EXPECT_TRUE(imStepIn.get() == 0 or
		imStepIn.get() == 1 or
		imStepIn.get() == 2);

    // the -1 is because step starts from 1
    const int index = (stepIn.get()-1) * 3 + imStepIn.get();
    EXPECT_TRUE( rhsIn.isApprox(rhs_.at(index)) );

    if (imStepIn.get() == 0){
      const IndepVarType expectedTime = (stepIn.get()-1)*dt_;
      EXPECT_TRUE( std::abs(expectedTime - timeIn) < 1e-13 );
    }
    if (imStepIn.get() == 1){
      const IndepVarType expectedTime = (stepIn.get()-1)*dt_ + dt_;
      EXPECT_TRUE( std::abs(expectedTime - timeIn) < 1e-13 );
    }
    if (imStepIn.get() == 2){
      const IndepVarType expectedTime = (stepIn.get()-1)*dt_ + dt_*0.5;
      EXPECT_TRUE( std::abs(expectedTime - timeIn) < 1e-13 );
    }
  }
};

#define COMMON_TEST(NAME, METHOD_SWITCH)			\
  using namespace pressio;					\
  std::map<int, Eigen::VectorXd> rhs;				\
  for (int i=0; i< 20; ++i){					\
    rhs[i] = Eigen::VectorXd::Random(3);			\
  }								\
  MyApp appObj(rhs);						\
  Eigen::VectorXd y(3);						\
  y.setZero();							\
  auto stepper = ode::create_##NAME##_stepper(appObj);		\
  StateObserver ob1;						\
  const double dt = 2.0;					\
  RhsObserver<METHOD_SWITCH> ob2(dt, rhs);			\
  auto bp = pressio::ode::steps_fixed_dt(0.0, pressio::ode::StepCount(2), dt); \
  ode::advance(stepper, y, bp, ob1, ob2);  \

TEST(ode_explicit_steppers, forward_euler){
  COMMON_TEST(forward_euler, 0);
}

TEST(ode_explicit_steppers, rk4){
  COMMON_TEST(rk4, 1);
}

TEST(ode_explicit_steppers, ab2){
  COMMON_TEST(ab2, 2);
}

TEST(ode_explicit_steppers, ssprk3){
  COMMON_TEST(ssprk3, 3);
}
