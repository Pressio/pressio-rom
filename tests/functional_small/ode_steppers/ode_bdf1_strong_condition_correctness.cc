
#include <gtest/gtest.h>
#include "pressio/ode_steppers.hpp"
#include <iomanip>

struct MyApp
{
  using independent_variable_type   = double;
  using state_type    = Eigen::VectorXd;
  using rhs_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

public:
  state_type createState() const{ return state_type(3); }

  rhs_type createRhs() const{ rhs_type f(3); return f; }

  jacobian_type createJacobian() const{
    jacobian_type J(3,3);
    // ensure that the diagonal elements exist
    for (int i=0; i<J.innerSize(); i++) {
      J.coeffRef(i, i) = 0;
    }
    return J;
  }

  void rhsAndJacobian(const state_type & /*unused*/,
		      const independent_variable_type& evaltime,
		      rhs_type & f,
		      std::optional<jacobian_type*> /*J*/) const{}
};

struct MyFakeSolver
{
  int count_={};
  int mustFailCount_={};

  MyFakeSolver(int v) : mustFailCount_(v){}

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    ++count_;
    std::cout << "SOLVE count = "  << count_ << std::endl;

    state(0) += 0.1;
    state(1) += 0.2;
    state(2) += 0.3;

    if (count_ == mustFailCount_){
      throw pressio::eh::NonlinearSolveFailure{};
    }
  }
};

// check that it we do one step of bdf1 and the solver fails,
// the solution is restored to what it was BEFORE attempting the solve

TEST(ode, implicit_bdf1_step_strong_condition_correctness_A)
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  using namespace pressio;
  using problem_t = MyApp;
  problem_t problemObj;
  auto stepper = ode::create_implicit_stepper(ode::StepScheme::BDF1,
					      problemObj);

  typename problem_t::state_type y(3); ops::fill(y, 1.);
  MyFakeSolver nls( 1 /*make it fail after 1 call */);

  // try to take one step, but since the solver is hardwired to fail
  // the ode stepper should also throw exception and we should
  // verify that the step y has not changed
  try{
    stepper(y, ode::StepStartAt<double>{0},
	    ode::StepCount{1},
	    ode::StepSize<double>{1.2},
	    nls);
  }
  catch(...){
    EXPECT_DOUBLE_EQ(y(0), 1.);
    EXPECT_DOUBLE_EQ(y(1), 1.);
    EXPECT_DOUBLE_EQ(y(2), 1.);
    std::cout << std::setprecision(14) << y << "\n";
  }

  PRESSIOLOG_FINALIZE();
}

TEST(ode, implicit_bdf1_step_strong_condition_correctness_B)
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  using namespace pressio;
  using problem_t = MyApp;
  problem_t problemObj;
  auto stepper = ode::create_implicit_stepper(ode::StepScheme::BDF1,
					      problemObj);

  typename problem_t::state_type y(3); ops::fill(y, 1.);

  MyFakeSolver nls( 3 /*make it fail after 3 calls */);

  // try to take one step, but since the solver is hardwired to fail
  // the ode stepper should also throw exception and we should
  // verify that the step y has not changed
  try{
    for (int i=0; i<3; ++i){
      stepper(y, ode::StepStartAt<double>{0},
	      ode::StepCount{i+1},
	      ode::StepSize<double>{1.2},
	      nls);
    }
  }
  catch(...){
    EXPECT_DOUBLE_EQ(y(0), 1.2);
    EXPECT_DOUBLE_EQ(y(1), 1.4);
    EXPECT_DOUBLE_EQ(y(2), 1.6);
    std::cout << std::setprecision(14) << y << "\n";
  }

  PRESSIOLOG_FINALIZE();
}
