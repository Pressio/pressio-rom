
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
		      std::optional<jacobian_type*> /*J*/) const
  {
    f.setConstant(1.);
  }
};

struct MyFakeSolver
{
  int count_={};
  int mustFailCount_={};
  Eigen::VectorXd resCopy_;

  MyFakeSolver(int v) : mustFailCount_(v), resCopy_(3){}

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    ++count_;
    std::cout << "SOLVE count = "  << count_ << std::endl;

    auto J = sys.createJacobian();
    sys.residualAndJacobian(state, resCopy_, std::optional<decltype(J)*>(&J));

    state(0) += 0.1;
    state(1) += 0.2;
    state(2) += 0.3;
    if (count_ == mustFailCount_){
      throw pressio::eh::NonlinearSolveFailure{};
    }
  }

  auto residualCopy() const{ return resCopy_; }
};

// check that it we do one step of bdf2 and the solver fails,
// the solution is restored to what it was BEFORE attempting the solve

TEST(ode, implicit_bdf2_step_strong_condition_correctness_A)
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  using namespace pressio;
  using problem_t = MyApp;
  problem_t problemObj;
  auto stepper = ode::create_implicit_stepper(ode::StepScheme::BDF2,
					      problemObj);

  typename problem_t::state_type y(3); ops::fill(y, 1.);
  MyFakeSolver nls( 1 /*make it fail after 1 call */);

  // try to take the first step, but since the solver is hardwired to fail
  // at the first call, the stepper should throw an exception and we should
  // verify that the state y has not changed.
  try{
    stepper(y, ode::StepStartAt<double>{0},
	    ode::StepCount{1}, ode::StepSize<double>{1.2},
	    nls);
  }
  catch(...){
    EXPECT_DOUBLE_EQ(y(0), 1.);
    EXPECT_DOUBLE_EQ(y(1), 1.);
    EXPECT_DOUBLE_EQ(y(2), 1.);
    std::cout << std::setprecision(14) << y << "\n";
  }

  // and after the failure is caught, if we try to repeat the same call
  stepper(y, ode::StepStartAt<double>{0},
	  ode::StepCount{1}, ode::StepSize<double>{1.2},
	  nls);

  // it is as if the step was never called so the residual computed should be:
  // R = - h * f
  auto r = nls.residualCopy();
  EXPECT_DOUBLE_EQ(r(0), -1.2);
  EXPECT_DOUBLE_EQ(r(1), -1.2);
  EXPECT_DOUBLE_EQ(r(2), -1.2);

  PRESSIOLOG_FINALIZE();
}

TEST(ode, implicit_bdf2_step_strong_condition_correctness_B)
{

  /*
    |------|------|------|
    y0    y1     y2     y3
    1     1.1
    1     1.2
    1     1.3
   */

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  using namespace pressio;
  using problem_t = MyApp;
  problem_t problemObj;
  auto stepper = ode::create_implicit_stepper(ode::StepScheme::BDF2,
					      problemObj);
  typename problem_t::state_type y(3); ops::fill(y, 1.);

  MyFakeSolver nls( 2 /*make it fail after these many calls */);

  try{
    for (int i=1; i<=2; ++i){
      stepper(y, ode::StepStartAt<double>{0},
	      ode::StepCount{i+1},
	      ode::StepSize<double>{1.2},
	      nls);
    }
  }
  catch(...){
    EXPECT_DOUBLE_EQ(y(0), 1.1);
    EXPECT_DOUBLE_EQ(y(1), 1.2);
    EXPECT_DOUBLE_EQ(y(2), 1.3);
    std::cout << std::setprecision(14) << y << "\n";
  }

  // and after the failure is caught, if we try to repeat the same call
  stepper(y, ode::StepStartAt<double>{0},
	  ode::StepCount{2},
	  ode::StepSize<double>{1.2},
	  nls);

  // it is as if the step was never called so the residual computed should be:
  // R = [1.1,1.2,1.3] -(4/3)*[1.1,1.2,1.3] -(1/3)*[1.,1.,1.] - h * f * (2/3)
  auto r = nls.residualCopy();
  EXPECT_DOUBLE_EQ(r(0), 1.1 - (4./3.)*1.1 + (1/3.) -1.2*(2/3.));
  EXPECT_DOUBLE_EQ(r(1), 1.2 - (4./3.)*1.2 + (1/3.) -1.2*(2/3.));
  EXPECT_DOUBLE_EQ(r(2), 1.3 - (4./3.)*1.3 + (1/3.) -1.2*(2/3.));
  std::cout << r << std::endl;

  PRESSIOLOG_FINALIZE();
}

TEST(ode, implicit_bdf2_step_strong_condition_correctness_C)
{

  /*
      Exception is thrown during this step
                      |
                      |
                      v
        1      2      3
    |------|------|------|
    y0    y1     y2     y3
    1     1.1   1.2     1.3
    1     1.2   1.4     1.6
    1     1.3   1.6     1.8
   */

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  using namespace pressio;
  using problem_t = MyApp;
  problem_t problemObj;
  auto stepper = ode::create_implicit_stepper(ode::StepScheme::BDF2,
					      problemObj);
  typename problem_t::state_type y(3); ops::fill(y, 1.);

  MyFakeSolver nls( 3 /*make it fail after these many calls */);
  try{
    for (int i=1; i<=3; ++i){
      stepper(y, ode::StepStartAt<double>{0},
	      ode::StepCount{i}, ode::StepSize<double>{1.2}, nls);
    }
  }
  catch(...){
    EXPECT_DOUBLE_EQ(y(0), 1.2);
    EXPECT_DOUBLE_EQ(y(1), 1.4);
    EXPECT_DOUBLE_EQ(y(2), 1.6);
    std::cout << std::setprecision(14) << y << "\n";
  }

  // and after the failure is caught, if we try to repeat the same call
  stepper(y, ode::StepStartAt<double>{0},
	  ode::StepCount{3}, ode::StepSize<double>{1.2}, nls);
  // it is as if the step was never called so the residual computed should be:
  // R = [1.2,1.4,1.6] -(4/3)*[1.2,1.4,1.6] -(1/3)*[1.1,1.2,1.3] - h * f * (2/3)
  auto r = nls.residualCopy();
  EXPECT_DOUBLE_EQ(r(0), 1.2 - (4./3.)*1.2 + (1/3.)*1.1 -1.2*(2/3.));
  EXPECT_DOUBLE_EQ(r(1), 1.4 - (4./3.)*1.4 + (1/3.)*1.2 -1.2*(2/3.));
  EXPECT_DOUBLE_EQ(r(2), 1.6 - (4./3.)*1.6 + (1/3.)*1.3 -1.2*(2/3.));
  std::cout << r << std::endl;

  PRESSIOLOG_FINALIZE();
}
