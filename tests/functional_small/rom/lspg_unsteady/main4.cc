
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_lspg_unsteady.hpp"

namespace{

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
  void solve(const SystemType & system, StateType & romState)
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
	system.residualAndJacobian(romState, R, &J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto ynp1 = ::pressio::ops::clone(R);
	compute_y(ynp1, romState);
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

	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	system.residualAndJacobian(romState, R, &J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto ynp1 = ::pressio::ops::clone(R);
	compute_y(ynp1, romState);

	StateType romState_n(3);
	romState_n(0) = 0; romState_n(1) = 1; romState_n(2) = 2;
	auto yn = ::pressio::ops::clone(R);
	compute_y(yn, romState_n);

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

	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
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
	system.residualAndJacobian(romState, R, &J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto ynp1 = ::pressio::ops::clone(R);
	compute_y(ynp1, romState);
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

	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	system.residualAndJacobian(romState, R, &J);
	// std::cout << "R = \n" << R << std::endl;
	// std::cout << "J = \n" << J << std::endl;

	auto ynp1 = ::pressio::ops::clone(R);
	compute_y(ynp1, romState);

	StateType romState_n(3);
	romState_n(0) = 2; romState_n(1) = 3; romState_n(2) = 4.;
	auto yn = ::pressio::ops::clone(R);
	compute_y(yn, romState_n);

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

	for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
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

class MyFom
{
  using phi_type = Eigen::MatrixXd;
  int N_ = {};

public:
  using time_type = double;
  using state_type = Eigen::VectorXd;
  using discrete_residual_type = state_type;

  MyFom(int N) : N_(N){}

  discrete_residual_type createDiscreteTimeResidual() const{
    return discrete_residual_type(N_);
  }

  phi_type createResultOfDiscreteTimeJacobianActionOn(const phi_type & B) const{
    return phi_type(N_, B.cols());
  }

  template<class StepCountType>
  void discreteTimeResidualAndJacobianAction(StepCountType,
					     double time,
					     double dt,
					     discrete_residual_type & R,
					     const phi_type & B,
               std::optional<phi_type*> JA,
					     const state_type & y_np1,
					     const state_type & y_n ) const
  {
    discrete_residual_type f(R.size());
    f.setZero();
    for (auto i=0; i<f.rows(); ++i){
     f(i) = y_np1(i) + time;
    }

    R = y_np1 -y_n - dt*f;

    if (bool(JA)){
      auto appJac = B;
      appJac.array() += time;
      *JA.value() = (B - dt*appJac);
    }
  }
};
}

TEST(rom_lspg_unsteady, test4)
{
  /* default lspg eigen with fully discrete API */

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

  auto problem = rom::lspg::create_unsteady_problem<2>(space, fomSystem);
  auto & stepper = problem; //lspgStepper();

  const double dt = 2.;
  FakeNonLinSolver<phi_t> nonLinSolver(N, phi, dt);
  Observer obs;
  auto policy = pressio::ode::steps_fixed_dt(0., pressio::ode::StepCount(2), dt);
  ode::advance(stepper, romState, policy, nonLinSolver, obs);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 4.);
  EXPECT_DOUBLE_EQ(romState[1], 5.);
  EXPECT_DOUBLE_EQ(romState[2], 6.);

  PRESSIOLOG_FINALIZE();
}
