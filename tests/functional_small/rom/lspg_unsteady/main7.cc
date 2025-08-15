
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_lspg_unsteady.hpp"

namespace{

struct MyFom
{
  using time_type    = double;
  using state_type     = Eigen::VectorXd;
  using rhs_type  = state_type;
  int N_ = {};

  MyFom(int N): N_(N){}

  rhs_type createRhs() const{
    auto R = rhs_type(N_);
    R.setConstant(0.);
    return R;
  }

  template<class OperandType>
  OperandType createResultOfJacobianActionOn(const OperandType & B) const
  {
    OperandType A(N_, B.cols());
    A.setConstant(0.);
    return A;
  }

  void rhs(const state_type & u,
	   time_type timeIn,
	   rhs_type & f) const
  {
    for (int i=0; i<f.rows(); ++i){
      f(i) = u(i) + timeIn;
    }
  }

  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     time_type time,
                     OperandType & A) const
  {
    A = B;
    A.array() += time;
  }
};

template<class phi_t>
struct FakeNonLinSolver
{
  int call_count_ = 0;
  int N_ = {};
  phi_t phi_;
  double dt_;

  FakeNonLinSolver(int N, phi_t phi, double dt)
    : N_(N), phi_(phi), dt_(dt){}

  template<class FomY_t, class RomStateType>
  void compute_y_fom(FomY_t & y, const RomStateType & yrom)
  {
    for (int i=0; i<N_; ++i){
      y(i) = 0;
      for (int j=0; j<3; ++j){
	y(i) += phi_(i,j) * yrom(j);
      }
    }
  }

  template<class FomY_t, class F_t>
  void compute_f(const FomY_t & y, F_t & f, double time)
  {
    for (int i=0; i<N_; ++i){
      f(i) = y(i) + time;
    }
  }

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    ++call_count_;
    auto R = system.createResidual();
    auto J = system.createJacobian();
    //using Jo_t = std::optional<decltype(J) *>;
    EXPECT_TRUE((std::size_t)pressio::ops::extent(R,0)==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,0)==(std::size_t)N_);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,1)==(std::size_t)3);
    std::cout << " CALL COUNT = " << call_count_ << '\n';

    //*******************************
    // call_count == 1
    // for this case, we should be doing BDF1
    // since this is the first step
    //*******************************
    if(call_count_==1){
      double predictionTime = 2.;
      //-----------------------
      // do solver iterator 1
      //-----------------------
      {
	EXPECT_NEAR(state(0), 0., 1e-15);
	EXPECT_NEAR(state(1), 1., 1e-15);
	EXPECT_NEAR(state(2), 2., 1e-15);

	system.residualAndJacobian(state, R, &J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto ynp1 = ::pressio::ops::clone(R);
	compute_y_fom(ynp1, state);
	auto f = ::pressio::ops::clone(R);
	compute_f(ynp1, f, predictionTime);
	// for callcount1, iter1: ynp1 and yn should be same,
	// for checking residual we can just use f only
	for (int i=0; i<N_; ++i){
	  EXPECT_NEAR(R(i), -dt_*f(i), 1e-12);
	}
	// and jacobian
	int count = 0;
	for (int i=0; i<N_; ++i){
	  for (int j=0; j<3; ++j){
	    const auto phi_val =  (double)count++;
	    EXPECT_NEAR(J(i,j), phi_val - dt_*(phi_val+predictionTime), 1e-12);
	  }
	}

	for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	EXPECT_NEAR(state(0), 1., 1e-15);
	EXPECT_NEAR(state(1), 2., 1e-15);
	EXPECT_NEAR(state(2), 3., 1e-15);

	system.residualAndJacobian(state, R, &J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto ynp1 = ::pressio::ops::clone(R);
	compute_y_fom(ynp1, state);
	auto f = ::pressio::ops::clone(R);
	compute_f(ynp1, f, predictionTime);

	StateType state_n(3);
	state_n(0) = 0; state_n(1) = 1; state_n(2) = 2;
	auto yn = ::pressio::ops::clone(R);
	compute_y_fom(yn, state_n);

	for (int i=0; i<N_; ++i){
	  EXPECT_NEAR(R(i), ynp1(i)-yn(i)-dt_*f(i), 1e-12);
	}
	int count = 0;
	for (int i=0; i<N_; ++i){
	  for (int j=0; j<3; ++j){
	    const auto phi_val =  (double)count++;
	    EXPECT_NEAR(J(i,j), phi_val - dt_*(phi_val+2.), 1e-12);
	  }
	}

	for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
      }
    }// callcount==1


    //*******************************
    // call_count == 2
    // for this case, we should be doing BDF2
    //*******************************
    if(call_count_==2){
      double predictionTime = 4.;
      //-----------------------
      // do solver iterator 1
      //-----------------------
      {
	EXPECT_NEAR(state(0), 2., 1e-15);
	EXPECT_NEAR(state(1), 3., 1e-15);
	EXPECT_NEAR(state(2), 4., 1e-15);

	system.residualAndJacobian(state, R, &J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto y_fom_np1 = ::pressio::ops::clone(R);
	compute_y_fom(y_fom_np1, state);

	StateType state_nm1(3);
	state_nm1 << 0., 1., 2.;
	auto y_fom_nm1 = ::pressio::ops::clone(R);
	compute_y_fom(y_fom_nm1, state_nm1);

	auto f = ::pressio::ops::clone(R);
	compute_f(y_fom_np1, f, predictionTime);

	for (int i=0; i<N_; ++i){
	  EXPECT_NEAR(R(i),
		      -(1./3.)*y_fom_np1(i) + (1./3.)*y_fom_nm1(i) -dt_*(2./3.)*f(i), 1e-13);
	}
	// and jacobian
	int count = 0;
	for (int i=0; i<N_; ++i){
	  for (int j=0; j<3; ++j){
	    const auto phi_val =  (double)count++;
	    EXPECT_NEAR(J(i,j), phi_val - (2./3.)*dt_*(phi_val+predictionTime), 1e-13);
	  }
	}

	for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	EXPECT_NEAR(state(0), 3., 1e-15);
	EXPECT_NEAR(state(1), 4., 1e-15);
	EXPECT_NEAR(state(2), 5., 1e-15);

	system.residualAndJacobian(state, R, &J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto y_fom_np1 = ::pressio::ops::clone(R);
	compute_y_fom(y_fom_np1, state);

	StateType state_n(3);
	state_n << 2., 3., 4.;
	auto y_fom_n = ::pressio::ops::clone(R);
	compute_y_fom(y_fom_n, state_n);

	StateType state_nm1(3);
	state_nm1 << 0., 1., 2.;
	auto y_fom_nm1 = ::pressio::ops::clone(R);
	compute_y_fom(y_fom_nm1, state_nm1);

	auto f = ::pressio::ops::clone(R);
	compute_f(y_fom_np1, f, predictionTime);

	for (int i=0; i<N_; ++i){
	  EXPECT_NEAR(R(i),
		      y_fom_np1(i) -(4./3.)*y_fom_n(i) + (1./3.)*y_fom_nm1(i) -dt_*(2./3.)*f(i), 1e-13);
	}
	int count = 0;
	for (int i=0; i<N_; ++i){
	  for (int j=0; j<3; ++j){
	    const auto phi_val =  (double)count++;
	    EXPECT_NEAR(J(i,j), phi_val - (2./3.)*dt_*(phi_val+predictionTime), 1e-12);
	  }
	}

	for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
      }
    }//callcount ==2


    //*******************************
    // call_count == 3
    // for this case, we should be doing BDF2
    //*******************************
    if(call_count_==3){
      double predictionTime = 6.;
      //-----------------------
      // do solver iterator 1
      //-----------------------
      {
	EXPECT_NEAR(state(0), 4., 1e-15);
	EXPECT_NEAR(state(1), 5., 1e-15);
	EXPECT_NEAR(state(2), 6., 1e-15);

	system.residualAndJacobian(state, R, &J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto y_fom_np1 = ::pressio::ops::clone(R);
	compute_y_fom(y_fom_np1, state);

	StateType state_nm1(3);
	state_nm1 << 2., 3., 4.;
	auto y_fom_nm1 = ::pressio::ops::clone(R);
	compute_y_fom(y_fom_nm1, state_nm1);

	auto f = ::pressio::ops::clone(R);
	compute_f(y_fom_np1, f, predictionTime);

	for (int i=0; i<N_; ++i){
	  EXPECT_NEAR(R(i),
		      -(1./3.)*y_fom_np1(i) + (1./3.)*y_fom_nm1(i) -dt_*(2./3.)*f(i), 1e-13);
	}
	// and jacobian
	int count = 0;
	for (int i=0; i<N_; ++i){
	  for (int j=0; j<3; ++j){
	    const auto phi_val =  (double)count++;
	    EXPECT_NEAR(J(i,j), phi_val - (2./3.)*dt_*(phi_val+predictionTime), 1e-13);
	  }
	}

	for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	EXPECT_NEAR(state(0), 5., 1e-15);
	EXPECT_NEAR(state(1), 6., 1e-15);
	EXPECT_NEAR(state(2), 7., 1e-15);

	system.residualAndJacobian(state, R, &J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto y_fom_np1 = ::pressio::ops::clone(R);
	compute_y_fom(y_fom_np1, state);

	StateType state_n(3);
	state_n << 4., 5., 6.;
	auto y_fom_n = ::pressio::ops::clone(R);
	compute_y_fom(y_fom_n, state_n);

	StateType state_nm1(3);
	state_nm1 << 2., 3., 4.;
	auto y_fom_nm1 = ::pressio::ops::clone(R);
	compute_y_fom(y_fom_nm1, state_nm1);

	auto f = ::pressio::ops::clone(R);
	compute_f(y_fom_np1, f, predictionTime);

	for (int i=0; i<N_; ++i){
	  EXPECT_NEAR(R(i),
		      y_fom_np1(i) -(4./3.)*y_fom_n(i) + (1./3.)*y_fom_nm1(i) -dt_*(2./3.)*f(i), 1e-13);
	}
	int count = 0;
	for (int i=0; i<N_; ++i){
	  for (int j=0; j<3; ++j){
	    const auto phi_val =  (double)count++;
	    EXPECT_NEAR(J(i,j), phi_val - (2./3.)*dt_*(phi_val+predictionTime), 1e-12);
	  }
	}

	for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
      }
    }//callcount ==3
  }
};

struct Observer{
  void operator()(pressio::ode::StepCount stepIn,
		  double /*time*/,
		  const Eigen::VectorXd & state) const
  {

    const auto step = stepIn.get();
    EXPECT_TRUE(step<=3);
    if (step==0){
      EXPECT_DOUBLE_EQ(state[0], 0.);
      EXPECT_DOUBLE_EQ(state[1], 1.);
      EXPECT_DOUBLE_EQ(state[2], 2.);
    }
    if (step==1){
      EXPECT_DOUBLE_EQ(state[0], 2.);
      EXPECT_DOUBLE_EQ(state[1], 3.);
      EXPECT_DOUBLE_EQ(state[2], 4.);
    }
    if (step==2){
      EXPECT_DOUBLE_EQ(state[0], 4.);
      EXPECT_DOUBLE_EQ(state[1], 5.);
      EXPECT_DOUBLE_EQ(state[2], 6.);
    }
    if (step==3){
      EXPECT_DOUBLE_EQ(state[0], 6.);
      EXPECT_DOUBLE_EQ(state[1], 7.);
      EXPECT_DOUBLE_EQ(state[2], 8.);
    }
  }
};
}

TEST(rom_lspg_unsteady, test7)
{
  /* default lspg eigen */

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  constexpr int N = 8;
  using fom_t = MyFom;
  fom_t fomSystem(N);

  using phi_t = Eigen::Matrix<double, -1,-1>;
  phi_t phi(N, 3);
  int count = 0;
  for (int i=0; i<N; ++i){
    for (int j=0; j<3; ++j){
      phi(i,j) = (double) count++;
    }
  }
  std::cout << phi << "\n";

  using namespace pressio;

  using reduced_state_type = Eigen::VectorXd;
  typename fom_t::state_type dummyFomState(N);
  constexpr bool isAffine = false;
  auto space = rom::create_trial_column_subspace<
    reduced_state_type>(phi, dummyFomState, isAffine);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  auto problem = rom::lspg::create_unsteady_problem(ode::StepScheme::BDF2,
						    space, fomSystem);
  auto & stepper = problem;

  const double dt = 2.;
  FakeNonLinSolver<phi_t> nonLinSolver(N, phi, dt);
  Observer obs;
  auto policy = pressio::ode::steps_fixed_dt(0., pressio::ode::StepCount(3), dt); 
  ode::advance(stepper, romState, policy, nonLinSolver, obs);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 6.);
  EXPECT_DOUBLE_EQ(romState[1], 7.);
  EXPECT_DOUBLE_EQ(romState[2], 8.);

  PRESSIOLOG_FINALIZE();
}
