
#include <gtest/gtest.h>
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"

struct MyApp2WithMM
{
  using independent_variable_type = double;
  using state_type           = Eigen::VectorXd;
  using rhs_type = state_type;
  using mass_matrix_type     = Eigen::MatrixXd;

  mutable int count1 = 0;
  mutable int count2 = 0;
  std::map<int, Eigen::VectorXd> & rhs_;
  std::map<int, Eigen::MatrixXd> & matrices_;

  MyApp2WithMM(std::map<int, Eigen::VectorXd> & rhs,
	       std::map<int, Eigen::MatrixXd> & matrices)
    : rhs_(rhs),
      matrices_(matrices){}

  state_type createState() const{
    state_type ret(3); ret.setZero();
    return ret;
  }

  rhs_type createRhs() const{
    rhs_type ret(3); ret.setZero();
    return ret;
  };

  mass_matrix_type createMassMatrix() const{
    mass_matrix_type ret(3,3); ret.setZero();
    return ret;
  };

  void massMatrixAndRhs(const state_type & /*unused*/,
			independent_variable_type evaltime,
			mass_matrix_type & M,
			rhs_type & rhs) const
  {
    rhs = Eigen::VectorXd::Random(rhs.rows());
    for (int i=0; i<rhs.size(); ++i){
      rhs(i) = evaltime;
    }

    M = Eigen::MatrixXd::Random(M.rows(), M.cols());
    for (int i=0; i<M.rows(); ++i){
      M(i,i) = evaltime;
    }

    rhs_[++count1] = rhs;
    matrices_[++count2] = M;
  };

};

struct MyApp2NoMM
{
  using independent_variable_type = double;
  using state_type           = Eigen::VectorXd;
  using rhs_type = state_type;

  mutable int count1 = 0;
  mutable int count2 = 0;
  const std::map<int, Eigen::VectorXd> & rhs_;
  const std::map<int, Eigen::MatrixXd> & matrices_;

  MyApp2NoMM(const std::map<int, Eigen::VectorXd> & rhs,
	     const std::map<int, Eigen::MatrixXd> & matrices)
    : rhs_(rhs), matrices_(matrices){}

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
    rhs = rhs_.at(++count1);
    // we want to make sure we have the correct time
    for (int i=0; i<rhs.size(); ++i){
      EXPECT_NEAR(rhs(i), evaltime, 1e-15);
    }

    auto M = createMassMatrix();
    massMatrix(y, evaltime, M);
    rhs = M.inverse()*rhs;
  };

private:
  using my_mass_matrix_type = Eigen::MatrixXd;

  my_mass_matrix_type createMassMatrix() const{
    my_mass_matrix_type ret(3,3); ret.setZero();
    return ret;
  };

  void massMatrix(const state_type & /*unused*/,
		  independent_variable_type evaltime,
		  my_mass_matrix_type & M) const
  {
    M = matrices_.at(++count2);
    for (int i=0; i<M.rows(); ++i){
      EXPECT_NEAR(M(i,i), evaltime, 1e-15);
    }
  };
};

struct LinearSolver2
{
  void solve(const Eigen::MatrixXd & A,
	     Eigen::VectorXd & x,
	     const Eigen::VectorXd & b)
  {
    x = A.colPivHouseholderQr().solve(b);
  }
};

#define ODE_ALL_SCHEME_MASS_MATRIX_CHECK_TEST(NAME)	\
  /* for given method, M dy/dt = f with mass-matrix API
     should be same as doing dy/dt = M^-1 f */			\
  std::cout << "\n";						\
  using namespace pressio;					\
  srand(342556331);						\
  const auto nsteps = ::pressio::ode::StepCount(4);		\
  const double dt = 2.;						\
  std::map<int, Eigen::VectorXd> rhs;				\
  std::map<int, Eigen::MatrixXd> massMatrices;			\
  /* first, solve problem using mass matrix API */		\
  Eigen::VectorXd y0(3);					\
  y0.setZero();{						\
    MyApp2WithMM appObj(rhs, massMatrices);			\
    auto stepperObj = ode::create_##NAME##_stepper(appObj);	\
    LinearSolver2 solver;					\
    y0(0) = 1.; y0(1) = 2.; y0(2) = 3.;					\
    auto bp = pressio::ode::steps_fixed_dt(0.0, nsteps, dt); \
    ode::advance(stepperObj, y0, bp, solver);	\
  }									\
  std::cout << "y0 : \n" << y0 << " \n";				\
  /* solve problem computing modified rhs using mass matrix inverse */	\
  Eigen::VectorXd y1(3);						\
  y1.setZero();{							\
    MyApp2NoMM appObj(rhs, massMatrices);				\
    auto stepperObj = ode::create_##NAME##_stepper(appObj);		\
    y1(0) = 1.; y1(1) = 2.; y1(2) = 3.;					\
    auto bp = pressio::ode::steps_fixed_dt(0.0, nsteps, dt); \
    ode::advance(stepperObj, y1, bp); \
  }									\
  std::cout << "y1 : \n" << y1 << " \n";				\
  EXPECT_NEAR( y0(0), y1(0), 1e-12);					\
  EXPECT_NEAR( y0(1), y1(1), 1e-12);					\
  EXPECT_NEAR( y0(2), y1(2), 1e-12);					\


TEST(ode_explicit_steppers, forward_euler_with_mass_matrix_use_inverse){
  ODE_ALL_SCHEME_MASS_MATRIX_CHECK_TEST(forward_euler)
}

TEST(ode_explicit_steppers, ab2_with_mass_matrix_use_inverse){
  ODE_ALL_SCHEME_MASS_MATRIX_CHECK_TEST(ab2)
}

TEST(ode_explicit_steppers, rk4_with_mass_matrix_use_inverse){
  ODE_ALL_SCHEME_MASS_MATRIX_CHECK_TEST(rk4)
}

TEST(ode_explicit_steppers, spprk3_with_mass_matrix_use_inverse){
  ODE_ALL_SCHEME_MASS_MATRIX_CHECK_TEST(ssprk3)
}
