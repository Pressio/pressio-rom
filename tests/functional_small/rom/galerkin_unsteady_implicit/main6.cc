
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_galerkin_unsteady.hpp"

namespace{

constexpr int N = 5;

using FomStateType = Eigen::VectorXd;
using FomRhsType = Eigen::VectorXd;
using RomStateType = Eigen::VectorXd;
using RomRType = Eigen::VectorXd;
using RomJType = Eigen::MatrixXd;

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

  Eigen::MatrixXd createResultOfMassMatrixActionOn(const Eigen::MatrixXd & operand) const{
    return Eigen::MatrixXd(N, operand.cols());
  }

  void applyMassMatrix(const Eigen::VectorXd & stateIn,
		       const Eigen::MatrixXd & operand,
		       double evalTime,
		       Eigen::MatrixXd & result) const
  {
    Eigen::MatrixXd M(N, N);
    M.setConstant(1);
    M.diagonal() = stateIn;
    result = M * operand;
  }
};

struct NonLinSolver
{
  int count_ = 0;

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    count_++;
    auto R = system.createResidual();
    auto J = system.createJacobian();
    //using Jo_t = std::optional<decltype(J) *>;

    // mimic iteration 1
    system.residualAndJacobian(state, R, &J);

    if (count_ == 1)
    {
      Eigen::VectorXd goldR(3); goldR << 0.,-70.,-140.;
      EXPECT_TRUE(R.isApprox(goldR));

      Eigen::MatrixXd goldJ(3,3);
      goldJ << 0.,0.,0., 20., 75., 130., 40., 150., 260.;
      EXPECT_TRUE(J.isApprox(goldJ));
    }
    for (int i=0; i<state.size(); ++i){ state[i] += 1.; }

    // mimic iteration 2
    system.residualAndJacobian(state, R, &J);
    if (count_ == 1)
    {
      Eigen::VectorXd goldR(3); goldR << 0.,80.,160.;
      EXPECT_TRUE(R.isApprox(goldR));
    }
    for (int i=0; i<state.size(); ++i){ state[i] += 1.; }

  }
};
}

TEST(rom_galerkin_implicit, default_with_massmatrix_bdf1)
{
  // implicit default galerkin with mass matrix

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
  NonLinSolver solver;

  auto policy = pressio::ode::steps_fixed_dt(time_type{0}, pressio::ode::StepCount(1), dt);
  pressio::ode::advance(problem, romState, policy, solver);

  PRESSIOLOG_FINALIZE();
}
