
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
  const std::vector<int> indices_ = {};
  int nstencil_ = {};

  MyFom(std::vector<int> ind, int nstencil)
  : N_(ind.size()), indices_(ind), nstencil_(nstencil){}

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
    EXPECT_TRUE((std::size_t)u.size()!=(std::size_t)f.size());
    EXPECT_TRUE((std::size_t)f.size()==(std::size_t)N_);
    for (std::size_t i=0; i<indices_.size(); ++i){
     f(i) = u(indices_[i]) + timeIn;
    }
  }

  template<class OperandType>
  void applyJacobian(const state_type & state,
                     const OperandType & B,
                     time_type time,
                     OperandType & A) const
  {
    EXPECT_TRUE((std::size_t)state.size()!=(std::size_t)N_);
    EXPECT_TRUE((std::size_t)A.rows()==(std::size_t)N_);

    for (std::size_t i=0; i<indices_.size(); ++i){
      for (int j=0; j< A.cols(); ++j){
        A(i,j) = B(indices_[i], j) + time;
      }
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

template<class ScalarType>
struct HypRedUpdaterEigen
{
  using vec_operand_type = Eigen::VectorXd;
  using mat_operand_type = Eigen::MatrixXd;
  const std::vector<int> rows_ = {};

  HypRedUpdaterEigen(std::vector<int> rows) : rows_(rows){}

  // a = alpha*a + beta*b (a,b potentially non with same distribution)
  void updateSampleMeshOperandWithStencilMeshOne(vec_operand_type & a, ScalarType alpha,
						 const vec_operand_type & b, ScalarType beta) const
  {
    for (std::size_t i=0; i<rows_.size(); ++i){
      a(i) = alpha*a(i) + beta*b(rows_[i]);
    }
  }

  void updateSampleMeshOperandWithStencilMeshOne(mat_operand_type & a, ScalarType alpha,
						 const mat_operand_type & b, ScalarType beta) const
  {
    for (std::size_t i=0; i<rows_.size(); ++i){
      for (int j=0; j<b.cols(); ++j){
	a(i,j) = alpha*a(i,j) + beta*b(rows_[i],j);
      }
    }
  }
};

struct Scaler{
  void operator()(const Eigen::VectorXd & statein,
		  double evalTime,
		  Eigen::VectorXd & residual,
		  std::optional<Eigen::MatrixXd *> jacobian) const
  {
    std::cout << "SCALING\n";
  }
};
}

TEST(rom_lspg_unsteady, test6)
{
  /* hyper-reduced lspg eigen for bdf1 with a trivial scaler
   should be identical to main1.cc.
   note that this WILL need to be changed to a non-trivial scaler
  */

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  const int nstencil = 15;
  const std::vector<int> sample_indices = {0,2,4,6,8,10,12,14};
  const int nSample  = sample_indices.size();
  MyFom fomSystem(sample_indices, nstencil);

  using phi_t = Eigen::Matrix<double, -1,-1>;
  phi_t phi(nstencil, 3);
  fill_phi(phi);
  using namespace pressio;

  using reduced_state_type = Eigen::VectorXd;
  typename MyFom::state_type dummyFomState(nstencil);
  constexpr bool isAffine = false;
  auto space = rom::create_trial_column_subspace<
    reduced_state_type>(phi, dummyFomState, isAffine);
  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  namespace plspg = rom::lspg::experimental;
  Scaler p;
  HypRedUpdaterEigen<double> comb(sample_indices);
  auto problem = plspg::create_unsteady_problem(ode::StepScheme::BDF1,
						space, fomSystem, comb, p);
  auto & stepper = problem; //lspgStepper();

  const double dt = 2.;
  FakeNonLinSolver nonLinSolver(nSample);
  Observer obs;
  auto policy = pressio::ode::steps_fixed_dt(0., pressio::ode::StepCount(2), dt);  
  ode::advance(stepper, romState, policy, nonLinSolver, obs);
  std::cout << romState << std::endl;
  EXPECT_DOUBLE_EQ(romState[0], 4.);
  EXPECT_DOUBLE_EQ(romState[1], 5.);
  EXPECT_DOUBLE_EQ(romState[2], 6.);

  PRESSIOLOG_FINALIZE();
}
