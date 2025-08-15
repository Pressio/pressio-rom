
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
    EXPECT_TRUE((std::size_t)u.size()==(std::size_t)f.size());
    EXPECT_TRUE((std::size_t)u.size()==(std::size_t)N_);

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
  void compute_y(FomY_t & y, const RomStateType & yrom)
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
    //
    // call_count == 1
    //
    //*******************************
    if(call_count_==1)
    {
      double mytime = 2.;
      //-----------------------
      // do solver iterator 1
      //-----------------------
      {
	system.residualAndJacobian(state, R, &J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto ynp1 = ::pressio::ops::clone(R);
	compute_y(ynp1, state);
	auto f = ::pressio::ops::clone(R);
	compute_f(ynp1, f, mytime);
	// for callcount1, iter1: ynp1 and yn should be same, so use f only
	for (int i=0; i<N_; ++i){
	  EXPECT_NEAR(R(i), -dt_*f(i), 1e-15);
	}
	int count = 0;
	for (int i=0; i<N_; ++i){
	  for (int j=0; j<3; ++j){
	    const auto phi_val =  (double)count++;
	    EXPECT_NEAR(J(i,j), phi_val - dt_*(phi_val+mytime), 1e-15);
	  }
	}

	for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	system.residualAndJacobian(state, R, &J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto ynp1 = ::pressio::ops::clone(R);
	compute_y(ynp1, state);

	StateType state_n(3);
	state_n(0) = 0; state_n(1) = 1; state_n(2) = 2;
	auto yn = ::pressio::ops::clone(R);
	compute_y(yn, state_n);

	auto f = ::pressio::ops::clone(R);
	compute_f(ynp1, f, mytime);

	for (int i=0; i<N_; ++i){
	  EXPECT_NEAR(R(i), ynp1(i)-yn(i)-dt_*f(i), 1e-15);
	}
	int count = 0;
	for (int i=0; i<N_; ++i){
	  for (int j=0; j<3; ++j){
	    const auto phi_val =  (double)count++;
	    EXPECT_NEAR(J(i,j), phi_val - dt_*(phi_val+2.), 1e-15);
	  }
	}

	for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
      }
    }

    //*******************************
    //
    // call_count == 2
    //
    //*******************************
    if(call_count_==2)
    {
      double mytime = 4.;
      //-----------------------
      // do solver iterator 1
      //-----------------------
      {
	system.residualAndJacobian(state, R, &J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto ynp1 = ::pressio::ops::clone(R);
	compute_y(ynp1, state);
	auto f = ::pressio::ops::clone(R);
	compute_f(ynp1, f, mytime);
	// for callcount2, iter1: ynp1 and yn should be same, so use f only
	for (int i=0; i<N_; ++i){
	  EXPECT_NEAR(R(i), -dt_*f(i), 1e-15);
	}
	int count = 0;
	for (int i=0; i<N_; ++i){
	  for (int j=0; j<3; ++j){
	    const auto phi_val =  (double)count++;
	    EXPECT_NEAR(J(i,j), phi_val - dt_*(phi_val+mytime), 1e-15);
	  }
	}

	for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	system.residualAndJacobian(state, R, &J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto ynp1 = ::pressio::ops::clone(R);
	compute_y(ynp1, state);

	StateType state_n(3);
	state_n(0) = 2; state_n(1) = 3; state_n(2) = 4.;
	auto yn = ::pressio::ops::clone(R);
	compute_y(yn, state_n);

	auto f = ::pressio::ops::clone(R);
	compute_f(ynp1, f, mytime);

	for (int i=0; i<N_; ++i){
	  EXPECT_NEAR(R(i), ynp1(i)-yn(i)-dt_*f(i), 1e-15);
	}
	int count = 0;
	for (int i=0; i<N_; ++i){
	  for (int j=0; j<3; ++j){
	    const auto phi_val =  (double)count++;
	    EXPECT_NEAR(J(i,j), phi_val - dt_*(phi_val+mytime), 1e-15);
	  }
	}

	for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
      }
    }
  }
};

struct Observer
{
  void operator()(pressio::ode::StepCount stepIn,
		  double /*time*/,
		  const Eigen::VectorXd & state) const
  {
    const auto step = stepIn.get();
    EXPECT_TRUE(step<=2);

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
  }
};

struct Scaler{
  void operator()(const Eigen::VectorXd & statein,
		  double evalTime,
		  Eigen::VectorXd & residual,
		  std::optional<Eigen::MatrixXd *> jacobian) const
  {
    std::cout << "SCALING\n";
  }
};
}

TEST(rom_lspg_unsteady, test5)
{
  /* default lspg eigen for bdf1 with a trivial scaler,
   should be identical to main1.cc.
   note that this WILL need to be changed to a non-trivial scale
  */

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

  Scaler scaler;
  auto problem = rom::lspg::experimental::create_unsteady_problem
    (pressio::ode::StepScheme::BDF1, space, fomSystem, scaler);
  auto & stepper = problem;//.lspgStepper();

  const double dt = 2.;
  FakeNonLinSolver<phi_t> nonLinSolver(N, phi, dt);
  Observer obs;
  auto policy = pressio::ode::steps_fixed_dt(0., pressio::ode::StepCount(2), dt);
  pressio::ode::advance(stepper, romState, policy, nonLinSolver, obs);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 4.);
  EXPECT_DOUBLE_EQ(romState[1], 5.);
  EXPECT_DOUBLE_EQ(romState[2], 6.);

  PRESSIOLOG_FINALIZE();
}
