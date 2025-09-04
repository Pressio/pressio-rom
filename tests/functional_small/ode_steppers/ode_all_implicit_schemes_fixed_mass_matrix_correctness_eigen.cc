
#include <gtest/gtest.h>
#include "pressio/ode_steppers.hpp"
#include "pressio/ode_advancers.hpp"

struct MyApp2WithMM
{
  using independent_variable_type = double;
  using state_type           = Eigen::VectorXd;
  using rhs_type = state_type;
  using mass_matrix_type     = Eigen::MatrixXd;
  using jacobian_type        = Eigen::MatrixXd;

  mutable int count1 = 0;
  mutable int count3 = 0;
  const std::map<int, rhs_type> & rhs_;
  const std::map<int, jacobian_type> & jacobians_;
  const mass_matrix_type & MM_;

  MyApp2WithMM(const std::map<int, rhs_type> & rhs,
	       const std::map<int, jacobian_type> & jacobians,
	       const mass_matrix_type & MMIn)
    : rhs_(rhs), jacobians_(jacobians), MM_(MMIn){}

  state_type createState() const{
    state_type ret(3); ret.setZero();
    return ret;
  }

  rhs_type createRhs() const{
    rhs_type ret(3); ret.setZero();
    return ret;
  };

  jacobian_type createJacobian() const{
    jacobian_type JJ(3,3);
    return JJ;
  };

  mass_matrix_type createMassMatrix() const{
    mass_matrix_type ret(3,3); ret.setZero();
    return ret;
  };

  void massMatrixAndRhsAndJacobian(const state_type & /*unused*/,
				   independent_variable_type /*unused*/,
				   mass_matrix_type & M,
				   rhs_type & rhs,
				   std::optional<jacobian_type*> JJ) const
  {
    rhs = rhs_.at(count1++);
    M = MM_;
    if (JJ){
      *(JJ.value()) = jacobians_.at(count3++);
    }
  };

};

struct MyApp2NoMM
{
  using independent_variable_type = double;
  using state_type           = Eigen::VectorXd;
  using rhs_type = state_type;
  using jacobian_type        = Eigen::MatrixXd;

  mutable int count1 = 0;
  mutable int count3 = 0;
  const std::map<int, Eigen::VectorXd> & rhs_;
  const Eigen::MatrixXd & MM_;
  const std::map<int, jacobian_type> & jacobians_;

  MyApp2NoMM(const std::map<int, Eigen::VectorXd> & rhs,
	     const Eigen::MatrixXd & MMIn,
	     const std::map<int, jacobian_type> & jacobians)
    : rhs_(rhs), MM_(MMIn), jacobians_(jacobians){}

  state_type createState() const{
    state_type ret(3); ret.setZero();
    return ret;
  }

  rhs_type createRhs() const{
    rhs_type ret(3); ret.setZero();
    return ret;
  };

  jacobian_type createJacobian() const{
    jacobian_type JJ(3,3);
    return JJ;
  };

  void rhsAndJacobian(const state_type & y,
		      independent_variable_type evaltime,
		      rhs_type & rhs,
		      std::optional<jacobian_type*> JJ) const
  {

    auto tmprhs = rhs_.at(count1++);
    rhs = MM_.inverse()*tmprhs;
    if(JJ){
      auto tmpJ = jacobians_.at(count3++);
      *(JJ.value()) = MM_.inverse()*tmpJ;
    }
  }
};

struct FakeNonLinearSolver1{
  int numFakeSolverIterations_{};
  int count_ = 0;
  const Eigen::MatrixXd & MM_;
  std::vector<Eigen::VectorXd> & residuals_;
  std::vector<Eigen::MatrixXd> & jacobians_;

  FakeNonLinearSolver1(int numFakeSolverIterations,
		       const Eigen::MatrixXd & MMIn,
		       std::vector<Eigen::VectorXd> & residuals,
		       std::vector<Eigen::MatrixXd> & jacobians)
    : numFakeSolverIterations_(numFakeSolverIterations),
      MM_(MMIn),
      residuals_(residuals),
      jacobians_(jacobians){}

  int count() const { return count_; }

  template<class SystemType, class StateType>
  void solve(SystemType & S, StateType & y)
  {
    auto R = S.createResidual();
    auto J = S.createJacobian();
    for (int i=0; i<numFakeSolverIterations_; ++i){
      count_++;

      S.residualAndJacobian(y, R, std::optional<decltype(J)*>(&J));

      auto MMinv = MM_.inverse();
      Eigen::VectorXd tmp = MMinv*R;
      residuals_.push_back(tmp);
      Eigen::MatrixXd tmp2 = MMinv*J;
      jacobians_.push_back(tmp2);

      for (int j=0; j<y.size(); ++j) y[j] += 1.;
    }
  }
};

struct FakeNonLinearSolver2{
  int count_ = 0;
  int numFakeSolverIterations_{};
  const std::vector<Eigen::VectorXd> & residuals_;
  std::vector<Eigen::MatrixXd> & jacobians_;

  FakeNonLinearSolver2(int numFakeSolverIterations,
		       const std::vector<Eigen::VectorXd> & residuals,
		       std::vector<Eigen::MatrixXd> & jacobians)
    : numFakeSolverIterations_(numFakeSolverIterations),
      residuals_(residuals),
      jacobians_(jacobians){}

  int count() const { return count_; }

  template<class SystemType, class StateType>
  void solve(SystemType & S, StateType & y)
  {
    auto R = S.createResidual();
    auto J = S.createJacobian();
    for (int i=0; i<numFakeSolverIterations_; ++i){
      S.residualAndJacobian(y, R, std::optional<decltype(J)*>(&J));
      EXPECT_TRUE( R.isApprox(residuals_[count_]) );
      EXPECT_TRUE( J.isApprox(jacobians_[count_++]) );

      for (int j=0; j<y.size(); ++j) y[j] += 1.;
    }
  }
};

#define ODE_MASS_MATRIX_CHECK_TEST(NAME)				\
  std::cout << "\n";							\
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::info);			\
  using namespace pressio;						\
  srand(342556331);							\
  const auto nsteps = ::pressio::ode::StepCount(4);			\
  const double dt = 2.;							\
  const int numFakeSolverIterations = 2;				\
  const int N = 3;							\
  /* fix MM once */							\
  const Eigen::MatrixXd MM = Eigen::MatrixXd::Random(3, 3);		\
  /**/									\
  std::map<int, Eigen::VectorXd> rhs;					\
  std::map<int, Eigen::MatrixXd> rhsJacobians;				\
  /**/									\
  /* to be safe, need to store many random instances */			\
  /* so that we don't get errors when running things */			\
  /* becuase we try to access instances that are not there */		\
  for (int i=0; i< 20; ++i){						\
    rhs[i] = Eigen::VectorXd::Random(N);				\
    rhsJacobians[i] = Eigen::MatrixXd::Random(N, N);			\
  }									\
  std::vector<Eigen::VectorXd> odeSchemeResidualsToCompare;		\
  std::vector<Eigen::MatrixXd> odeSchemeJacobiansToCompare;		\
  /**/									\
  /* solve problem using mass matrix API */				\
  /**/									\
  Eigen::VectorXd y0(3);						\
  y0.setZero();								\
  {									\
    MyApp2WithMM appObj(rhs, rhsJacobians, MM);				\
    auto stepperObj = ode::create_##NAME##_stepper(appObj);		\
    y0(0) = 1.; y0(1) = 2.; y0(2) = 3.;					\
    FakeNonLinearSolver1 solver(numFakeSolverIterations, MM,		\
				odeSchemeResidualsToCompare,		\
				odeSchemeJacobiansToCompare);		\
    auto policy = ode::steps_fixed_dt(0.0, nsteps, dt); \
    ode::advance(stepperObj, y0, policy, solver);	\
    EXPECT_TRUE(solver.count()== nsteps.get()*numFakeSolverIterations); \
  }									\
  std::cout << "y0 : \n" << y0 << " \n";				\
  /**/									\
  /* solver problem using normal API */					\
  /**/									\
  Eigen::VectorXd y1(3);						\
  y1.setZero();								\
  {									\
    MyApp2NoMM appObj(rhs, MM, rhsJacobians);				\
    auto stepperObj = ode::create_##NAME##_stepper(appObj);		\
    y1(0) = 1.; y1(1) = 2.; y1(2) = 3.;					\
    FakeNonLinearSolver2 solver(numFakeSolverIterations,		\
				odeSchemeResidualsToCompare,		\
				odeSchemeJacobiansToCompare);		\
    auto policy = ode::steps_fixed_dt(0.0, nsteps, dt); \
    ode::advance(stepperObj, y1, policy, solver);	\
    EXPECT_TRUE(solver.count() == nsteps.get()*numFakeSolverIterations); \
  }									\
  std::cout << "y1 : \n" << y1 << " \n";				\
  /**/									\
  EXPECT_NEAR( y0(0), y1(0), 1e-12);					\
  EXPECT_NEAR( y0(1), y1(1), 1e-12);					\
  EXPECT_NEAR( y0(2), y1(2), 1e-12);					\
  PRESSIOLOG_FINALIZE();						\

TEST(ode_implicit_steppers, bdf1_with_fixed_mass_matrix_use_inverse)
{
  /*
    How can we check if, e.g., BDF1 with MM works correctly?

    Recall that for BDF1 with MM the residual is:
       R = M_n+1( y_n+1 - y_n) -h*f_n+1

    Suppose that we knew M_n+1^-1 at each step, then
    using the MM API we compute at each step
      g1 = M_n+1^-1 R
         = ( y_n+1 - y_n) -h*M_n+1^-1*f_n+1

    which should match what we get bu solving a problem
    WITHOUT MASS MATRIX but with a "modified rhs",
    where we compute:
      R2 = ( y_n+1 - y_n) -h*M_n+1^-1*f_n+1

    The test uses fake (random) data for the operators.
    And we use a fake solver so that we can check that the
    operators match, i.e. g1 == R2

    So in other other words we exercise the two distinct APIs,
    i.e. regular and mass matrix APIs, for implicit steppers,
    in such a way that we know we should get matching results.
    Note that this is can only be done for
    some schemes, not everything.
   */

  ODE_MASS_MATRIX_CHECK_TEST(bdf1)
}

TEST(ode_implicit_steppers, bdf2_with_fixed_mass_matrix_use_inverse){
  // similar explanation as for BDF1
  ODE_MASS_MATRIX_CHECK_TEST(bdf2)
}

// TEST(ode_implicit_steppers, cn_with_fixed_mass_matrix_use_inverse){
//   ODE_MASS_MATRIX_CHECK_TEST(cranknicolson)
// }
