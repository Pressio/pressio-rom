
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_galerkin_unsteady.hpp"

namespace{

constexpr int N = 7;

using FomStateType = Eigen::VectorXd;
using FomRhsType = Eigen::VectorXd;
using RomStateType = Eigen::VectorXd;
using RomRType = Eigen::VectorXd;
using RomJType = Eigen::MatrixXd;

struct Gold
{
  RomRType R_Step1i1;
  RomJType J_Step1i1;
  RomRType R_Step1i2;
  RomJType J_Step1i2;

  Gold() :
    R_Step1i1(3), J_Step1i1(3,3),
    R_Step1i2(3), J_Step1i2(3,3)
  {
    R_Step1i1 << 0., -98., -196;
    R_Step1i2 << 1., 1.-140., 1.-280.;

    J_Step1i1.row(0) << 1.,   0.,  0.;
    J_Step1i1.row(1) << 28., 43., 56.;
    J_Step1i1.row(2) << 56., 84., 113.;
    J_Step1i2 = J_Step1i1;
  }
};

struct MyFom
{
  using time_type = double;
  using state_type = FomStateType;
  using rhs_type = FomRhsType;

  MyFom(){}

  rhs_type createRhs() const{
    rhs_type r(N);
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

  Eigen::MatrixXd createResultOfJacobianActionOn(const Eigen::MatrixXd & A) const{
    return Eigen::MatrixXd(N, A.cols());
  }

  void applyJacobian(const state_type & s,
		     const Eigen::MatrixXd & A,
		     const time_type & evaltime,
		     Eigen::MatrixXd & result) const
  {
    result.col(0).setConstant(-2.);
    result.col(1).setConstant(-3.);
    result.col(2).setConstant(-4.);
  }
};

struct NonLinearSolver
{
  Gold gold;
  int stepTracker = 0;

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    stepTracker++;

    auto R = system.createResidual();
    auto J = system.createJacobian();
    //using Jo_t = std::optional<decltype(J) *>;

    //
    // do fake first iteration
    //
    system.residualAndJacobian(state, R, &J);
    std::cout << R << "\n";
    std::cout << J << "\n";
    if (stepTracker == 1){
      EXPECT_TRUE( R.isApprox(gold.R_Step1i1) );
      EXPECT_TRUE( J.isApprox(gold.J_Step1i1) );
    }
    // we fake an update to the state
    for (int i=0; i<state.size(); ++i){ state[i] += 1.; }

    //
    // do fake second iteration
    //
    system.residualAndJacobian(state, R, &J);
    std::cout << R << "\n";
    if (stepTracker == 1){
      EXPECT_TRUE( R.isApprox(gold.R_Step1i2) );
      EXPECT_TRUE( J.isApprox(gold.J_Step1i2) );
    }
    // we fake an update to the state
    for (int i=0; i<state.size(); ++i){ state[i] += 1.; }
  }
};
}

TEST(rom_galerkin_implicit, default_bdf1)
{
  // test for default implicit galerkin using BDF1
  // all numbers have been computed manually

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  // create fom
  using fom_t = MyFom;
  fom_t fomSystem;

  // create trial space
  using basis_t = Eigen::MatrixXd;
  basis_t phi(N, 3);
  phi.col(0).setConstant(0.);
  phi.col(1).setConstant(1.);
  phi.col(2).setConstant(2.);

  typename MyFom::state_type shift(N);
  auto space = pressio::rom::create_trial_column_subspace<RomStateType>(phi,
								 shift,
								 false);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  const auto odeScheme = pressio::ode::StepScheme::BDF1;
  namespace gal = pressio::rom::galerkin;
  auto problem = gal::create_unsteady_implicit_problem(odeScheme, space, fomSystem);

  using time_type = typename fom_t::time_type;
  const time_type dt = 2.;
  NonLinearSolver solver;

  auto policy = pressio::ode::steps_fixed_dt(time_type{0}, pressio::ode::StepCount(1), dt);
  pressio::ode::advance(problem, romState, policy, solver);
  std::cout << romState << std::endl;
  EXPECT_TRUE(romState[0] == 2.);
  EXPECT_TRUE(romState[1] == 3.);
  EXPECT_TRUE(romState[2] == 4.);

  PRESSIOLOG_FINALIZE();
}
