
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
  using time_type       = double;
  using state_type        = Eigen::VectorXd;
  using rhs_type     = state_type;

  std::vector<int> sampleMeshIndices_;
  int nstencil_ = {};
  int nsample_ = {};

  MyFom(int nstencil, const std::vector<int> & sampleMeshIndices)
    : sampleMeshIndices_(sampleMeshIndices),
      nstencil_(nstencil),
      nsample_(sampleMeshIndices.size()){}

  rhs_type createRhs() const{
    return rhs_type(nsample_);
  }

  void rhs(const state_type & u,
	   const time_type timeIn,
	   rhs_type & f) const
  {
    EXPECT_TRUE((std::size_t)u.size()!=(std::size_t)f.size());
    EXPECT_TRUE((std::size_t)f.size()==(std::size_t) nsample_);

    for (std::size_t i=0; i<sampleMeshIndices_.size(); ++i){
     f(i) = u(sampleMeshIndices_[i]) + timeIn;
    }
  }
};

class HypRedOperator
{
  using operator_type = Eigen::MatrixXd;
  operator_type matrix_;

public:
  using time_type = double;
  using rhs_operand_type = Eigen::VectorXd;

  HypRedOperator(const operator_type & phi) : matrix_(phi){}

  void operator()(const rhs_operand_type & operand,
		  time_type timein,
		  Eigen::VectorXd & result) const
  {
    result = matrix_.transpose() * operand;
  }
};
}

TEST(rom_galerkin_explicit, hyperreduced_velo_euler_forward)
{
  /*
    hypred velocity Galerkin, Euler forward:
    this is prtty much the same as main1 except that
    we fake hypereduction by mimicing a stencil/sample mesh
  */

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  const std::vector<int> sampleMeshIndices = {1,3,5,7,9,11,13,15,17,19};
  const int nstencil = 20;

  MyFom fomSystem(nstencil, sampleMeshIndices);

  using basis_t = Eigen::MatrixXd;
  basis_t phi(nstencil, 3);
  phi.col(0).setConstant(0.);
  phi.col(1).setConstant(1.);
  phi.col(2).setConstant(2.);
  // to check that we are doing things right, we
  // corrupt all locations that are NOT part of the sample mesh
  // so that if things are right the numbers are not used
  phi.row(0).setConstant(-111.);
  phi.row(2).setConstant(-111.);
  phi.row(4).setConstant(111.);
  phi.row(6).setConstant(423.);
  phi.row(8).setConstant(-21.);
  phi.row(10).setConstant(423.);
  phi.row(12).setConstant(-21.);
  phi.row(14).setConstant(423.);
  phi.row(16).setConstant(-21.);
  phi.row(18).setConstant(-21.);

  using reduced_state_type = Eigen::VectorXd;
  typename MyFom::state_type shift(nstencil);
  auto space = pressio::rom::create_trial_column_subspace<reduced_state_type>(phi, shift, false);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  /* hypred operator must be applicable to an operand on *sample* mesh */
  basis_t matForProj(sampleMeshIndices.size(), 3);
  matForProj.col(0).setConstant(0.);
  matForProj.col(1).setConstant(1.);
  matForProj.col(2).setConstant(2.);
  HypRedOperator hrOp(matForProj);

  const auto odeScheme = pressio::ode::StepScheme::ForwardEuler;
  namespace gal = pressio::rom::galerkin;
  auto problem = gal::create_unsteady_explicit_problem(odeScheme, space, fomSystem, hrOp);

  const double dt = 1.;
  Observer obs;
  auto policy = pressio::ode::steps_fixed_dt(0., pressio::ode::StepCount(2), dt);
  pressio::ode::advance(problem, romState, policy, obs);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 0.);
  EXPECT_DOUBLE_EQ(romState[1], 2611.);
  EXPECT_DOUBLE_EQ(romState[2], 5222.);

  PRESSIOLOG_FINALIZE();
}
