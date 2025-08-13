
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
  const std::vector<int> indices_to_corrupt_ = {};

  MyFom(int N, std::vector<int> ind) : N_(N), indices_to_corrupt_(ind){}

  rhs_type createRhs() const{ return rhs_type(N_); }

  template<class OperandType>
  OperandType createResultOfJacobianActionOn(const OperandType & B) const
  {
    OperandType A(N_, B.cols());
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
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     f(it) = -1114;
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
    // corrupt some to ensure masking works
    for (auto & it : indices_to_corrupt_){
     A.row(it).setConstant(-1114);
    }
  }
};

struct FakeNonLinSolver
{
  int call_count_ = 0;
  int N_ = {};

  FakeNonLinSolver(int N) : N_(N){}

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

    //*******************************
    //
    // call_count == 1
    //
    //*******************************
    if(call_count_==1)
    {
      //-----------------------
      // do solver iterator 1
      //-----------------------
      {
	system.residualAndJacobian(state, R, &J);
	//std::cout << "S " << call_count_ << " \n" << R << std::endl;
	//std::cout << "S " << call_count_ << " \n" << J << std::endl;

	EXPECT_DOUBLE_EQ(R(0), -14.);
	EXPECT_DOUBLE_EQ(R(1), -32.);
	EXPECT_DOUBLE_EQ(R(2), -50.);
	EXPECT_DOUBLE_EQ(R(3), -68.);
	EXPECT_DOUBLE_EQ(R(4), -86.);
	EXPECT_DOUBLE_EQ(R(5), -104.);
	EXPECT_DOUBLE_EQ(R(6), -122.);
	EXPECT_DOUBLE_EQ(R(7), -140.);

	EXPECT_DOUBLE_EQ(J(0,0), -4.);
	EXPECT_DOUBLE_EQ(J(0,1), -5.);
	EXPECT_DOUBLE_EQ(J(0,2), -6.);
	EXPECT_DOUBLE_EQ(J(1,0), -7.);
	EXPECT_DOUBLE_EQ(J(1,1), -8.);
	EXPECT_DOUBLE_EQ(J(1,2), -9.);
	EXPECT_DOUBLE_EQ(J(2,0), -10.);
	EXPECT_DOUBLE_EQ(J(2,1), -11.);
	EXPECT_DOUBLE_EQ(J(2,2), -12.);
	EXPECT_DOUBLE_EQ(J(3,0), -13.);
	EXPECT_DOUBLE_EQ(J(3,1), -14.);
	EXPECT_DOUBLE_EQ(J(3,2), -15.);
	EXPECT_DOUBLE_EQ(J(4,0), -16.);
	EXPECT_DOUBLE_EQ(J(4,1), -17.);
	EXPECT_DOUBLE_EQ(J(4,2), -18.);
	EXPECT_DOUBLE_EQ(J(5,0), -19.);
	EXPECT_DOUBLE_EQ(J(5,1), -20.);
	EXPECT_DOUBLE_EQ(J(5,2), -21.);
	EXPECT_DOUBLE_EQ(J(6,0), -22.);
	EXPECT_DOUBLE_EQ(J(6,1), -23.);
	EXPECT_DOUBLE_EQ(J(6,2), -24.);
	EXPECT_DOUBLE_EQ(J(7,0), -25.);
	EXPECT_DOUBLE_EQ(J(7,1), -26.);
	EXPECT_DOUBLE_EQ(J(7,2), -27.);

	for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	system.residualAndJacobian(state, R, &J);
        //std::cout << "S " << call_count_ << " \n" << R << std::endl;
	//std::cout << "S " << call_count_ << " \n" << J << std::endl;

	EXPECT_DOUBLE_EQ(R(0), -17.);
	EXPECT_DOUBLE_EQ(R(1), -44.);
	EXPECT_DOUBLE_EQ(R(2), -71.);
	EXPECT_DOUBLE_EQ(R(3), -98.);
	EXPECT_DOUBLE_EQ(R(4), -125.);
	EXPECT_DOUBLE_EQ(R(5), -152.);
	EXPECT_DOUBLE_EQ(R(6), -179.);
	EXPECT_DOUBLE_EQ(R(7), -206.);

	EXPECT_DOUBLE_EQ(J(0,0), -4.);
	EXPECT_DOUBLE_EQ(J(0,1), -5.);
	EXPECT_DOUBLE_EQ(J(0,2), -6.);
	EXPECT_DOUBLE_EQ(J(1,0), -7.);
	EXPECT_DOUBLE_EQ(J(1,1), -8.);
	EXPECT_DOUBLE_EQ(J(1,2), -9.);
	EXPECT_DOUBLE_EQ(J(2,0), -10.);
	EXPECT_DOUBLE_EQ(J(2,1), -11.);
	EXPECT_DOUBLE_EQ(J(2,2), -12.);
	EXPECT_DOUBLE_EQ(J(3,0), -13.);
	EXPECT_DOUBLE_EQ(J(3,1), -14.);
	EXPECT_DOUBLE_EQ(J(3,2), -15.);
	EXPECT_DOUBLE_EQ(J(4,0), -16.);
	EXPECT_DOUBLE_EQ(J(4,1), -17.);
	EXPECT_DOUBLE_EQ(J(4,2), -18.);
	EXPECT_DOUBLE_EQ(J(5,0), -19.);
	EXPECT_DOUBLE_EQ(J(5,1), -20.);
	EXPECT_DOUBLE_EQ(J(5,2), -21.);
	EXPECT_DOUBLE_EQ(J(6,0), -22.);
	EXPECT_DOUBLE_EQ(J(6,1), -23.);
	EXPECT_DOUBLE_EQ(J(6,2), -24.);
	EXPECT_DOUBLE_EQ(J(7,0), -25.);
	EXPECT_DOUBLE_EQ(J(7,1), -26.);
	EXPECT_DOUBLE_EQ(J(7,2), -27.);

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
      //-----------------------
      // do solver iterator 1
      //-----------------------
      {
	system.residualAndJacobian(state, R, &J);
	//std::cout << "S " << call_count_ << " \n" << R << std::endl;
	//std::cout << "S " << call_count_ << " \n" << J << std::endl;

	EXPECT_DOUBLE_EQ(R(0), -30.);
	EXPECT_DOUBLE_EQ(R(1), -84.);
	EXPECT_DOUBLE_EQ(R(2), -138.);
	EXPECT_DOUBLE_EQ(R(3), -192.);
	EXPECT_DOUBLE_EQ(R(4), -246.);
	EXPECT_DOUBLE_EQ(R(5), -300.);
	EXPECT_DOUBLE_EQ(R(6), -354.);
	EXPECT_DOUBLE_EQ(R(7), -408.);

	EXPECT_DOUBLE_EQ(J(0,0), -8.);
	EXPECT_DOUBLE_EQ(J(0,1), -9.);
	EXPECT_DOUBLE_EQ(J(0,2), -10.);
	EXPECT_DOUBLE_EQ(J(1,0), -11.);
	EXPECT_DOUBLE_EQ(J(1,1), -12.);
	EXPECT_DOUBLE_EQ(J(1,2), -13.);
	EXPECT_DOUBLE_EQ(J(2,0), -14.);
	EXPECT_DOUBLE_EQ(J(2,1), -15.);
	EXPECT_DOUBLE_EQ(J(2,2), -16.);
	EXPECT_DOUBLE_EQ(J(3,0), -17.);
	EXPECT_DOUBLE_EQ(J(3,1), -18.);
	EXPECT_DOUBLE_EQ(J(3,2), -19.);
	EXPECT_DOUBLE_EQ(J(4,0), -20.);
	EXPECT_DOUBLE_EQ(J(4,1), -21.);
	EXPECT_DOUBLE_EQ(J(4,2), -22.);
	EXPECT_DOUBLE_EQ(J(5,0), -23.);
	EXPECT_DOUBLE_EQ(J(5,1), -24.);
	EXPECT_DOUBLE_EQ(J(5,2), -25.);
	EXPECT_DOUBLE_EQ(J(6,0), -26.);
	EXPECT_DOUBLE_EQ(J(6,1), -27.);
	EXPECT_DOUBLE_EQ(J(6,2), -28.);
	EXPECT_DOUBLE_EQ(J(7,0), -29.);
	EXPECT_DOUBLE_EQ(J(7,1), -30.);
	EXPECT_DOUBLE_EQ(J(7,2), -31.);

	for (int i=0; i<state.size(); ++i){ state(i) += 1.; }
      }

      //-----------------------
      // do solver iterator 2
      //-----------------------
      {
	system.residualAndJacobian(state, R, &J);
        //std::cout << "S " << call_count_ << " \n" << R << std::endl;
	//std::cout << "S " << call_count_ << " \n" << J << std::endl;

	EXPECT_DOUBLE_EQ(R(0), -33.);
	EXPECT_DOUBLE_EQ(R(1), -96.);
	EXPECT_DOUBLE_EQ(R(2), -159.);
	EXPECT_DOUBLE_EQ(R(3), -222.);
	EXPECT_DOUBLE_EQ(R(4), -285.);
	EXPECT_DOUBLE_EQ(R(5), -348.);
	EXPECT_DOUBLE_EQ(R(6), -411.);
	EXPECT_DOUBLE_EQ(R(7), -474.);

	EXPECT_DOUBLE_EQ(J(0,0), -8.);
	EXPECT_DOUBLE_EQ(J(0,1), -9.);
	EXPECT_DOUBLE_EQ(J(0,2), -10.);
	EXPECT_DOUBLE_EQ(J(1,0), -11.);
	EXPECT_DOUBLE_EQ(J(1,1), -12.);
	EXPECT_DOUBLE_EQ(J(1,2), -13.);
	EXPECT_DOUBLE_EQ(J(2,0), -14.);
	EXPECT_DOUBLE_EQ(J(2,1), -15.);
	EXPECT_DOUBLE_EQ(J(2,2), -16.);
	EXPECT_DOUBLE_EQ(J(3,0), -17.);
	EXPECT_DOUBLE_EQ(J(3,1), -18.);
	EXPECT_DOUBLE_EQ(J(3,2), -19.);
	EXPECT_DOUBLE_EQ(J(4,0), -20.);
	EXPECT_DOUBLE_EQ(J(4,1), -21.);
	EXPECT_DOUBLE_EQ(J(4,2), -22.);
	EXPECT_DOUBLE_EQ(J(5,0), -23.);
	EXPECT_DOUBLE_EQ(J(5,1), -24.);
	EXPECT_DOUBLE_EQ(J(5,2), -25.);
	EXPECT_DOUBLE_EQ(J(6,0), -26.);
	EXPECT_DOUBLE_EQ(J(6,1), -27.);
	EXPECT_DOUBLE_EQ(J(6,2), -28.);
	EXPECT_DOUBLE_EQ(J(7,0), -29.);
	EXPECT_DOUBLE_EQ(J(7,1), -30.);
	EXPECT_DOUBLE_EQ(J(7,2), -31.);

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

template<class T>
void fill_phi(T & phi)
{
  /*
    - phi in R^{10,3}:
        phi[0,:]=0,1,2
        phi[1,:]=-1,-1,-1
        phi[2,:]=3,4,5
        phi[3,:]=-1,-1,-1
        phi[4,:]=6,7,8
        phi[5,:]=-1,-1,-1
        phi[6,:]=9,10,11
        phi[7,:]=-1,-1,-1
        phi[8,:]=12,13,14
        phi[9,:]=-1,-1,-1
        phi[10,:]=15,16,17
        phi[11,:]=-1,-1,-1
        phi[12,:]=18,19,20
        phi[13,:]=-1,-1,-1
        phi[14,:]=21,22,23
  */

  using sc_t = typename pressio::Traits<T>::scalar_type;
  const int nrows = pressio::ops::extent(phi, 0);
  const int ncols = pressio::ops::extent(phi, 1);
  int count = 0;
  for (int i=0; i<nrows; ++i){
    for (int j=0; j<ncols; ++j){
      if (i % 2 == 0){
        phi(i,j) = (sc_t) count++;
      }
      else{
        phi(i,j) = (sc_t) -1;
      }
    }
  }
}


class MyMasker
{
  const std::vector<int> sample_indices_ = {};

public:
  MyMasker(std::vector<int> sample_indices) : sample_indices_(sample_indices){}

  auto createResultOfMaskActionOn(const Eigen::VectorXd & /*operand*/) const{
    return Eigen::VectorXd(sample_indices_.size());
  }

  auto createResultOfMaskActionOn(const Eigen::MatrixXd & operand) const{
    return Eigen::MatrixXd(sample_indices_.size(), operand.cols());
  }

  template<class T1, class T2>
  void operator()(const Eigen::MatrixBase<T1> & operand,
		  Eigen::MatrixBase<T2> & result) const
  {
    for (std::size_t i=0; i<sample_indices_.size(); ++i){
      for (int j=0; j<operand.cols(); ++j){
        result(i,j) = operand(sample_indices_[i],j);
      }
    }
  }
};

}

TEST(rom_lspg_unsteady, test3)
{
  /* masked lspg eigen */

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  const std::vector<int> rows_to_corrupt_ = {1,3,5,7,9,11,13};\
  const std::vector<int> sample_indices = {0,2,4,6,8,10,12,14};\
  const int nMasked = sample_indices.size();\
  const int N = (int) (rows_to_corrupt_.size() + sample_indices.size());\
  MyFom fomSystem(N, rows_to_corrupt_);

  using phi_t = Eigen::Matrix<double, -1,-1>;
  phi_t phi(N, 3);
  fill_phi(phi);

  using namespace pressio;

  using reduced_state_type = Eigen::VectorXd;
  typename MyFom::state_type dummyFomState(N);
  constexpr bool isAffine = false;
  auto space = rom::create_trial_column_subspace<
    reduced_state_type>(phi, dummyFomState, isAffine);
  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  MyMasker masker(sample_indices);
  auto problem = rom::lspg::create_unsteady_problem
    (ode::StepScheme::BDF1, space, fomSystem, masker);
  auto & stepper = problem; //lspgStepper();

  const double dt = 2.;
  FakeNonLinSolver nonLinSolver(nMasked);
  Observer obs;
  auto policy = pressio::ode::steps_fixed_dt(0., pressio::ode::StepCount(2), dt); 
  ode::advance(stepper, romState, policy, nonLinSolver, obs);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 4.);
  EXPECT_DOUBLE_EQ(romState[1], 5.);
  EXPECT_DOUBLE_EQ(romState[2], 6.);

  PRESSIOLOG_FINALIZE();
}
