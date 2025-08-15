
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_galerkin_unsteady.hpp"

namespace{

using phi_t = Eigen::Matrix<double, -1,-1>;
constexpr double dt = 2.;
constexpr int numModes = 3;

auto create_basis(int N){
  phi_t phi(N, numModes);
  int count = 0;
  for (int i=0; i<N; ++i){
    for (int j=0; j<numModes; ++j){
      phi(i,j) = (double) count++;
    }
  }
  return phi;
}

template<class phi_t>
struct FakeNonLinSolver
{
  int call_count_ = 0;
  int N_ = {};
  phi_t phi_;
  double dt_;

  FakeNonLinSolver(int N, phi_t phi, double dt)
    : N_(N), phi_(phi), dt_(dt){}

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    const auto phi= create_basis(N_);

    ++call_count_;
    std::cout << " CALL COUNT = " << call_count_ << '\n';

    auto R = system.createResidual();
    auto J = system.createJacobian();
    EXPECT_TRUE((std::size_t)pressio::ops::extent(R,0)==(std::size_t)3);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,0)==(std::size_t)3);
    EXPECT_TRUE((std::size_t)pressio::ops::extent(J,1)==(std::size_t)3);
    //using Jo_t = std::optional<decltype(J) *>;

    //-----------------------
    // fake a solver iterator 1
    //-----------------------
    system.residualAndJacobian(state, R, &J);

    if (call_count_ == 1){
      const double predictionTime = 2.;

      {
	// check residual
	Eigen::VectorXd ynp1(numModes); ynp1 << 0.,1.,2.;
	Eigen::VectorXd yn(numModes);   yn   << 0.,1.,2.;
	Eigen::VectorXd ynp1_fom = phi*ynp1;
	Eigen::VectorXd yn_fom   = phi*yn;
	Eigen::VectorXd f_fom    = ynp1_fom;
	f_fom.array() += predictionTime;
	Eigen::VectorXd R_fom = ynp1_fom - yn_fom - dt*f_fom;
	Eigen::VectorXd goldR = phi.transpose() * R_fom;
	EXPECT_TRUE(R.isApprox(goldR));
      }

      {
	// check jacobian
	Eigen::MatrixXd tmp = phi;
	tmp.array() += predictionTime;
	Eigen::MatrixXd fomJA = phi - dt*tmp;
	Eigen::MatrixXd goldJ = phi.transpose() * fomJA;
	EXPECT_TRUE(goldJ.isApprox(J));
      }

    }
    else if (call_count_ == 2)
    {
      const double predictionTime = 4.;

      {
	Eigen::VectorXd ynp1(numModes); ynp1 << 2.,3.,4.;
	Eigen::VectorXd yn(numModes);   yn   << 2.,3.,4.;
	Eigen::VectorXd ynp1_fom = phi*ynp1;
	Eigen::VectorXd yn_fom   = phi*yn;
	Eigen::VectorXd f_fom    = ynp1_fom;
	f_fom.array() += predictionTime;
	Eigen::VectorXd R_fom = ynp1_fom - yn_fom - dt*f_fom;
	Eigen::VectorXd goldR = phi.transpose() * R_fom;
	EXPECT_TRUE(R.isApprox(goldR));
      }

      {
	// check jacobian
	Eigen::MatrixXd tmp = phi;
	tmp.array() += predictionTime;
	Eigen::MatrixXd fomJA = phi - dt*tmp;
	Eigen::MatrixXd goldJ = phi.transpose() * fomJA;
	EXPECT_TRUE(goldJ.isApprox(J));
      }
    }

    for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

    //-----------------------
    // fake a solver iterator 2
    //-----------------------
    system.residualAndJacobian(state, R, &J);

    if (call_count_ == 1){
      const double predictionTime = 2.;

      {
	Eigen::VectorXd ynp1(numModes); ynp1 << 1.,2.,3.;
	Eigen::VectorXd yn(numModes);   yn   << 0.,1.,2.;
	Eigen::VectorXd ynp1_fom = phi*ynp1;
	Eigen::VectorXd yn_fom   = phi*yn;
	Eigen::VectorXd f_fom    = ynp1_fom;
	f_fom.array() += predictionTime;
	Eigen::VectorXd R_fom = ynp1_fom - yn_fom - dt*f_fom;
	Eigen::VectorXd goldR = phi.transpose() * R_fom;
	EXPECT_TRUE(R.isApprox(goldR));
      }

      {
	// check jacobian
	Eigen::MatrixXd tmp = phi;
	tmp.array() += predictionTime;
	Eigen::MatrixXd fomJA = phi - dt*tmp;
	Eigen::MatrixXd goldJ = phi.transpose() * fomJA;
	EXPECT_TRUE(goldJ.isApprox(J));
      }

    }
    else if (call_count_ == 2){
      const double predictionTime = 4.;

      {
	Eigen::VectorXd ynp1(numModes); ynp1 << 3.,4.,5.;
	Eigen::VectorXd yn(numModes);   yn   << 2.,3.,4.;
	Eigen::VectorXd ynp1_fom = phi*ynp1;
	Eigen::VectorXd yn_fom   = phi*yn;
	Eigen::VectorXd f_fom    = ynp1_fom;
	f_fom.array() += predictionTime;
	Eigen::VectorXd R_fom = ynp1_fom - yn_fom - dt*f_fom;
	Eigen::VectorXd goldR = phi.transpose() * R_fom;
	EXPECT_TRUE(R.isApprox(goldR));
      }

      {
	// check jacobian
	Eigen::MatrixXd tmp = phi;
	tmp.array() += predictionTime;
	Eigen::MatrixXd fomJA = phi - dt*tmp;
	Eigen::MatrixXd goldJ = phi.transpose() * fomJA;
	EXPECT_TRUE(goldJ.isApprox(J));
      }
    }

    for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
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
    EXPECT_TRUE(state.size()==numModes);

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

    R = y_np1 - y_n - dt*f;

    if (bool(JA)){
      auto appJac = B;
      appJac.array() += time;
      *JA.value() = (B - dt*appJac);
    }
  }
};
}

TEST(rom_galerkin_implicit, default_fullydiscrete_n2)
{
  /* default galerkin impliacit eigen with fully discrete API */

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  constexpr int N = 5;
  using fom_t = MyFom;
  fom_t fomSystem(N);

  auto phi = create_basis(N);
  std::cout << phi << "\n";

  using reduced_state_type = Eigen::VectorXd;
  typename fom_t::state_type dummyFomState(N);
  constexpr bool isAffine = false;
  auto space = pressio::rom::create_trial_column_subspace<
    reduced_state_type>(phi, dummyFomState, isAffine);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  auto problem = pressio::rom::galerkin::create_unsteady_implicit_problem<2>(space, fomSystem);

  FakeNonLinSolver<phi_t> nonLinSolver(N, phi, dt);
  Observer obs;
  auto policy = pressio::ode::steps_fixed_dt(0., pressio::ode::StepCount(2), dt);
  pressio::ode::advance(problem, romState, policy, nonLinSolver, obs);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 4.);
  EXPECT_DOUBLE_EQ(romState[1], 5.);
  EXPECT_DOUBLE_EQ(romState[2], 6.);

  PRESSIOLOG_FINALIZE();
}
