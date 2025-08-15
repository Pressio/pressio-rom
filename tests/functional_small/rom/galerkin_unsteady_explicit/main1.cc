
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_galerkin_unsteady.hpp"

namespace{
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
      EXPECT_DOUBLE_EQ(state[0], 0.);
      EXPECT_DOUBLE_EQ(state[1], 51.);
      EXPECT_DOUBLE_EQ(state[2], 102.);
    }
    if (step==2){
      EXPECT_DOUBLE_EQ(state[0], 0.);
      EXPECT_DOUBLE_EQ(state[1], 2611.);
      EXPECT_DOUBLE_EQ(state[2], 5222.);
    }
  }
};

struct MyFom
{
  using time_type = double;
  using state_type = Eigen::VectorXd;
  using rhs_type = state_type;
  int N_ = {};

  MyFom(int N): N_(N){}

  rhs_type createRhs() const{
    rhs_type r(N_);
    r.setConstant(0);
    return r;
  }

  void rhs(const state_type & u,
	   const time_type evalTime,
	   rhs_type & f) const
  {
    for (decltype(f.rows()) i=0; i<f.rows(); ++i){
      f(i) = u(i) + evalTime;
    }
  }
};
}

TEST(rom_galerkin_explicit, default)
{
  /*
    default Galerkin, Euler forward

    - run two steps: t0 -> t1 -> t2
    - dt = 1.

    - phi in R^{10,3}:
        phi[:,0]=0, phi[:,0]=1, phi[:,0]=2

    - initial romState = [0,1,2]

    - fom velocity f(y,t) = y+t always

    step0:
      rom state is the initial condition

    step1:
      y_fom = phi*[0,1,2]
      rom_state|_step1 = [0,1,2]^T + phi^T f(y_fom, t=0.) = [0,51,102]

    step2:
      y_fom = phi*[0,51,102] = [255, ..., 255]^T
      rom_state|_step2 = [0,51,102]^T + phi^T f(y_fom, t=1.) = [0, 2611, 5222]
  */

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  // create fom
  constexpr int N = 10;
  using fom_t = MyFom;
  fom_t fomSystem(N);

  // create trial space
  using basis_t = Eigen::MatrixXd;
  basis_t phi(N, 3);
  phi.col(0).setConstant(0.);
  phi.col(1).setConstant(1.);
  phi.col(2).setConstant(2.);

  using reduced_state_type = Eigen::VectorXd;
  typename fom_t::state_type shift(N);
  auto space = pressio::rom::create_trial_column_subspace<reduced_state_type>(phi, shift, false);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  const auto odeScheme = pressio::ode::StepScheme::ForwardEuler;
  namespace gal = pressio::rom::galerkin;
  auto problem = gal::create_unsteady_explicit_problem(odeScheme, space, fomSystem);

  using time_type = typename fom_t::time_type;
  const time_type dt = 1.;
  Observer obs;  
  auto policy = pressio::ode::steps_fixed_dt(time_type{0}, pressio::ode::StepCount(2), dt);  
  pressio::ode::advance(problem, romState, policy, obs);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 0.);
  EXPECT_DOUBLE_EQ(romState[1], 2611.);
  EXPECT_DOUBLE_EQ(romState[2], 5222.);

  PRESSIOLOG_FINALIZE();
}
