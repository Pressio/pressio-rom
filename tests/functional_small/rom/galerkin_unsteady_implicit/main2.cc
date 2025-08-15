
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_galerkin_unsteady.hpp"

namespace{

struct MyFom
{
  using time_type       = double;
  using state_type        = Eigen::VectorXd;
  using rhs_type = state_type;

  int N_ = {};
  int nStencil_ ={};
  const std::vector<int> indices = {1,3,5,7,9,11,13,15,17,19};

  MyFom(int N, int nStencil) : N_(N), nStencil_(nStencil){
    EXPECT_TRUE((std::size_t)N==(std::size_t)indices.size());
  }

  rhs_type createRhs() const{ return rhs_type(N_); }

  template<class OperandType>
  OperandType createResultOfJacobianActionOn(const OperandType & B) const{
    OperandType A(N_, B.cols());
    return A;
  }

  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     const time_type & timeIn,
                     OperandType & A) const
  {
    EXPECT_TRUE((std::size_t)state.size()!=(std::size_t)N_);
    EXPECT_TRUE((std::size_t)A.rows()==(std::size_t)N_);

    for (std::size_t i=0; i<indices.size(); ++i){
      for (int j=0; j< A.cols(); ++j){
        A(i,j) = B(indices[i], j);
      }
    }
    A.array() += timeIn;
  }

  void rhs(const state_type & u,
	   const time_type timeIn,
	   rhs_type & f) const
  {
    EXPECT_TRUE((std::size_t)u.size()!=(std::size_t)f.size());
    EXPECT_TRUE((std::size_t)f.size()==(std::size_t)N_);
    for (std::size_t i=0; i<indices.size(); ++i){
     f(i) = u(indices[i]) + timeIn;
    }
  }
};

class HypRedOperator
{
  using operator_type = Eigen::MatrixXd;
  operator_type matrix_;

public:
  using time_type = double;
  using rhs_operand_type = typename MyFom::rhs_type;
  using jacobian_action_operand_type = Eigen::MatrixXd;

  HypRedOperator(const operator_type & phi) : matrix_(phi){}

  void operator()(const rhs_operand_type & operand,
		  time_type timein,
		  Eigen::VectorXd & result) const
  {
    result = matrix_.transpose() * operand;
  }

  void operator()(const jacobian_action_operand_type & operand,
		  time_type timein,
		  Eigen::MatrixXd & result) const
  {
    result = matrix_.transpose() * operand;
  }
};

struct NonLinSolver
{
  int call_count_ = 0;

  template<class SystemType, class StateType>
  void solve(const SystemType & system, StateType & state)
  {
    ++call_count_;
    auto R = system.createResidual();
    auto J = system.createJacobian();
    //using Jo_t = std::optional<decltype(J) *>;

    //
    // call_count == 1
    //
    if(call_count_==1)
    {
      // do solver iterator 1
      system.residualAndJacobian(state, R, &J);
      // std::cout << "S " << call_count_ << " \n" << R << std::endl;
      // std::cout << "S " << call_count_ << " \n" << J << std::endl;
      EXPECT_DOUBLE_EQ(R[0], 0.);
      EXPECT_DOUBLE_EQ(R[1], -140.);
      EXPECT_DOUBLE_EQ(R[2], -280.);

      EXPECT_DOUBLE_EQ(J(0,0),   1.);
      EXPECT_DOUBLE_EQ(J(0,1),   0.);
      EXPECT_DOUBLE_EQ(J(0,2),   0.);
      EXPECT_DOUBLE_EQ(J(1,0), -40.);
      EXPECT_DOUBLE_EQ(J(1,1), -59.);
      EXPECT_DOUBLE_EQ(J(1,2), -80.);
      EXPECT_DOUBLE_EQ(J(2,0), -80.);
      EXPECT_DOUBLE_EQ(J(2,1),-120.);
      EXPECT_DOUBLE_EQ(J(2,2),-159.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

      // do solver iterator 2
      system.residualAndJacobian(state, R, &J);
      EXPECT_DOUBLE_EQ(R[0], 1.);
      EXPECT_DOUBLE_EQ(R[1], -199.);
      EXPECT_DOUBLE_EQ(R[2], -399.);

      EXPECT_DOUBLE_EQ(J(0,0),   1.);
      EXPECT_DOUBLE_EQ(J(0,1),   0.);
      EXPECT_DOUBLE_EQ(J(0,2),   0.);
      EXPECT_DOUBLE_EQ(J(1,0), -40.);
      EXPECT_DOUBLE_EQ(J(1,1), -59.);
      EXPECT_DOUBLE_EQ(J(1,2), -80.);
      EXPECT_DOUBLE_EQ(J(2,0), -80.);
      EXPECT_DOUBLE_EQ(J(2,1),-120.);
      EXPECT_DOUBLE_EQ(J(2,2),-159.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
    }

    //
    // call_count == 2
    //
    if(call_count_==2)
    {
      // do solver iterator 1
      system.residualAndJacobian(state, R, &J);
      EXPECT_DOUBLE_EQ(R[0],    0.);
      EXPECT_DOUBLE_EQ(R[1], -300.);
      EXPECT_DOUBLE_EQ(R[2], -600.);

      EXPECT_DOUBLE_EQ(J(0,0),   1.);
      EXPECT_DOUBLE_EQ(J(0,1),   0.);
      EXPECT_DOUBLE_EQ(J(0,2),   0.);
      EXPECT_DOUBLE_EQ(J(1,0), -80.);
      EXPECT_DOUBLE_EQ(J(1,1), -99.);
      EXPECT_DOUBLE_EQ(J(1,2), -120.);
      EXPECT_DOUBLE_EQ(J(2,0), -160.);
      EXPECT_DOUBLE_EQ(J(2,1), -200.);
      EXPECT_DOUBLE_EQ(J(2,2), -239.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }

      // do solver iterator 2
      system.residualAndJacobian(state, R, &J);
      EXPECT_DOUBLE_EQ(R[0],    1.);
      EXPECT_DOUBLE_EQ(R[1], -359.);
      EXPECT_DOUBLE_EQ(R[2], -719.);

      EXPECT_DOUBLE_EQ(J(0,0),   1.);
      EXPECT_DOUBLE_EQ(J(0,1),   0.);
      EXPECT_DOUBLE_EQ(J(0,2),   0.);
      EXPECT_DOUBLE_EQ(J(1,0), -80.);
      EXPECT_DOUBLE_EQ(J(1,1), -99.);
      EXPECT_DOUBLE_EQ(J(1,2), -120.);
      EXPECT_DOUBLE_EQ(J(2,0), -160.);
      EXPECT_DOUBLE_EQ(J(2,1), -200.);
      EXPECT_DOUBLE_EQ(J(2,2), -239.);

      for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
    }
  }
};
}

TEST(rom_galerkin_implicit, hyperreduced_bdf1)
{
  // test for hypred implicit galerkin using BDF1
  // all numbers have been computed manually

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  const int nStencil = 20;
  const int nSample  = 10;
  MyFom fomSystem(nSample, nStencil);
  typename MyFom::state_type fomReferenceState(nStencil);

  using phi_t = Eigen::MatrixXd;
  phi_t phi(nStencil, 3);
  phi.col(0).setConstant(0.);
  phi.col(1).setConstant(1.);
  phi.col(2).setConstant(2.);
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
  typename MyFom::state_type shift(nStencil);
  auto space = pressio::rom::create_trial_column_subspace<reduced_state_type>(phi,
								       shift,
								       false);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  phi_t matForProj(nSample, 3);
  matForProj.col(0).setConstant(0.);
  matForProj.col(1).setConstant(1.);
  matForProj.col(2).setConstant(2.);
  HypRedOperator hrOp(matForProj);

  const auto odeScheme = pressio::ode::StepScheme::BDF1;
  namespace gal = pressio::rom::galerkin;
  auto problem = gal::create_unsteady_implicit_problem(odeScheme, space, fomSystem, hrOp);

  const double dt = 2.;
  NonLinSolver solver;
  auto policy = pressio::ode::steps_fixed_dt(0., pressio::ode::StepCount(2), dt);
  pressio::ode::advance(problem, romState, policy, solver);
  std::cout << romState << std::endl;
  EXPECT_TRUE(romState[0] == 4.);
  EXPECT_TRUE(romState[1] == 5.);
  EXPECT_TRUE(romState[2] == 6.);

  PRESSIOLOG_FINALIZE();
}
