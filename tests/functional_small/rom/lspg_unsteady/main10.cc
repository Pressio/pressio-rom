
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_lspg_unsteady.hpp"

namespace{

template<class phi_t>
struct FakeNonLinSolver
{
  int N_ = {};
  phi_t phi_;
  double dt_;

  FakeNonLinSolver(int N, phi_t phi, double dt)
    : N_(N), phi_(phi), dt_(dt){}

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & romState)
  {
    auto R = system.createResidual();
    auto J = system.createJacobian();

    // ensure that we call the system multiples times to simulate multiple
    // iterations of a nonlinear solver
    for (int k=0; k<11; ++k){
      system.residualAndJacobian(romState, R, &J);
      for (int i=0; i<romState.size(); ++i){ romState(i) += 1.; }
    }
    // std::cout << "R = \n" << R << std::endl;
    // std::cout << "J = \n" << J << std::endl;
  }
};

class MyFom
{
  using phi_type = Eigen::MatrixXd;
  int N_ = {};
  mutable int hookCount_ = 0;

public:
  using time_type = double;
  using state_type = Eigen::VectorXd;
  using discrete_residual_type = state_type;

  explicit MyFom(int N) : N_(N){}

  discrete_residual_type createDiscreteTimeResidual() const{
    return discrete_residual_type(N_);
  }

  phi_type createResultOfDiscreteTimeJacobianActionOn(const phi_type & B) const{
    return phi_type(N_, B.cols());
  }

  template<class StepIntType>
  void preStepHook(StepIntType sc, double time, double dt,
               const state_type &,
               const state_type &) const
  {
    // this should be called only once per time step
    PRESSIOLOG_INFO("hook: {} {} {}\n", sc, time, dt);

    EXPECT_TRUE(sc <= 3);
    EXPECT_DOUBLE_EQ(dt, 2.);
    hookCount_++;
    if (hookCount_ == 1){
      EXPECT_TRUE(sc == 1);
      EXPECT_DOUBLE_EQ(time, 0.);
    }
    else if (hookCount_ == 2){
      EXPECT_TRUE(sc == 2);
      EXPECT_DOUBLE_EQ(time, 2.);
    }
    else if (hookCount_ == 3){
      EXPECT_TRUE(sc == 3);
      EXPECT_DOUBLE_EQ(time, 4.);
    }
  }

  template<class StepIntType>
  void discreteTimeResidualAndJacobianAction(StepIntType,
					     double time,
					     double dt,
					     discrete_residual_type & R,
					     const phi_type & B,
					     std::optional<phi_type*> JA,
					     const state_type & y_np1,
					     const state_type & y_n ) const
  {
    // no op since we don't really care about setting values here
  }
};
}

TEST(rom_lspg_unsteady, test10)
{
  /*
    default lspg eigen with fully discrete API
    where we have a callback hook that needs to be called
    before every new time step but only once per time step
  */

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  constexpr int N = 8;
  using fom_t = MyFom;
  fom_t fomSystem(N);

  using phi_t = Eigen::Matrix<double, -1,-1>;
  phi_t phi(N, 3);
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
  auto & stepper = problem.lspgStepper();

  const double dt = 2.;
  FakeNonLinSolver<phi_t> nonLinSolver(N, phi, dt);
  ode::advance_n_steps(stepper, romState, 0., dt,
		       ode::StepCount(3), nonLinSolver);
  std::cout << romState << std::endl;
  // EXPECT_DOUBLE_EQ(romState[0], 4.);
  // EXPECT_DOUBLE_EQ(romState[1], 5.);
  // EXPECT_DOUBLE_EQ(romState[2], 6.);

  PRESSIOLOG_FINALIZE();
}
