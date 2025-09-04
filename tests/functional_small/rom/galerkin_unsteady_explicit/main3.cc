
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_galerkin_unsteady.hpp"

namespace{

struct MyFom
{
  using time_type       = double;
  using state_type        = Eigen::VectorXd;
  using rhs_type     = state_type;
  int N_ = {};
  const std::vector<int> indices_to_corrupt_ = {};

  MyFom(int N, std::vector<int> ind)
    : N_(N), indices_to_corrupt_(ind){}

  rhs_type createRhs() const{
    return rhs_type(N_);
  }

  void rhs(const state_type & u,
	   const time_type timeIn,
	   rhs_type & f) const
  {
    for (auto i=0; i<f.rows(); ++i){
     f(i) = u(i) + timeIn;
    }
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     f(it) = -1114;
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

struct MyMasker
{
  const std::vector<int> sample_indices_ = {};
  using operand_type = typename MyFom::rhs_type;
  using result_type = Eigen::VectorXd;

  MyMasker(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  result_type createResultOfMaskActionOn(const operand_type & operand) const{
    return result_type(sample_indices_.size());
  }

  void operator()(const operand_type & operand, result_type & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      result(i) = operand(sample_indices_[i]);
    }
  }
};
}

TEST(rom_galerkin_explicit, masked_velo_euler_forward)
{
  /*
    check correctness of masked Galerkin explicit
    this test, *after the mask is applied*,
    is doing the same thing the default galerkin main1.cc
  */

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  constexpr int nFull = 20;\
  const std::vector<int> sample_indices = {0,2,4,6,8,10,12,14,16,18};
  const int nMasked = sample_indices.size();
  /* corrupt indices are those that we mess up on purpose */
  const std::vector<int> corrupt_indices = {1,7,13,19};
  MyFom fomSystem(nFull, corrupt_indices);

  using phi_t = Eigen::MatrixXd;
  phi_t phiFull(nFull, 3);
  phiFull.col(0).setConstant(0.);
  phiFull.col(1).setConstant(1.);
  phiFull.col(2).setConstant(2.);
  /* corrupting the entries in phiFull is only way we can ensure the masking works properly */
  phiFull.row(1).setConstant(-111.);
  phiFull.row(7).setConstant(111.);
  phiFull.row(11).setConstant(423.);
  phiFull.row(17).setConstant(-21.);

  using reduced_state_type = Eigen::VectorXd;
  typename MyFom::state_type shift(nFull);
  auto space = pressio::rom::create_trial_column_subspace<reduced_state_type>(phiFull, shift, false);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  /* hypred operator must be applicable to an operand on *sample* mesh */
  phi_t matForProj(nMasked, 3);
  for (int i = 0; i < nMasked; ++i){
    matForProj(i, 0) = phiFull(sample_indices[i],0);
    matForProj(i, 1) = phiFull(sample_indices[i],1);
    matForProj(i, 2) = phiFull(sample_indices[i],2);
  }
  HypRedOperator hrOp(matForProj);

  MyMasker masker(sample_indices);

  const auto odeScheme = pressio::ode::StepScheme::ForwardEuler;
  namespace gal = pressio::rom::galerkin;
  auto problem = gal::create_unsteady_explicit_problem(odeScheme, space,
						       fomSystem,
						       masker, hrOp);

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
