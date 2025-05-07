
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_lspg_unsteady.hpp"
#include <random>

namespace{

constexpr int N = 5;
using __this_test_phi_type = Eigen::MatrixXd;
using __this_test_rom_state_type = Eigen::VectorXd;
using __this_test_fom_state_type = Eigen::VectorXd;
constexpr double minStepSize = 0.1;
constexpr double reductionFactor = 2.;

bool __this_test_failureHappended = false;

// if the min step size = 0.1 and desired time horizon is 2,
// we need to have enough states stored
std::map<int, std::vector<double>> romStateAtStepStarting =
  {
    {0,  {0,0,0}},
    {1,  {0,1,2}},
    {2,  {5,6,7}},
    {3,  {10,11,12}},
    {4,  {15,16,17}},
    {5,  {20,21,22}},
    {6,  {25,26,27}},
    {7,  {30,31,32}},
    {8,  {35,36,37}},
    {9,  {40,41,42}},
    {10, {45,46,47}},
    {11, {50,51,52}},
    {12, {55,56,57}},
    {13, {60,61,62}},
    {14, {65,66,67}},
    {15, {70,71,72}},
    {16, {75,76,77}},
    {17, {80,81,82}},
    {18, {85,86,87}},
    {19, {90,91,92}},
    {20, {95,96,97}},
    {21, {100,101,102}},
    {22, {105,106,107}}
  };

auto create_gold_phi(){
  __this_test_phi_type phi(N,3);
  int count = 0;
  for (int i=0; i<N; ++i){
    for (int j=0; j<3; ++j){
      phi(i,j) = (double) count++;
    }
  }
  return phi;
}

auto stdvec_to_eigen(const std::vector<double> & v){
  Eigen::VectorXd a(v.size());
  for (int i=0; i<(int)v.size(); ++i){ a(i) = v[i]; }
  return a;
}

//
// FOM CLASS
//
class MyFom
{
public:
  using time_type = double;
  using state_type = Eigen::VectorXd;
  using discrete_residual_type = state_type;
  using dist_type = std::uniform_real_distribution<double>;

  mutable int count_ = 0;
  std::random_device rd;
  mutable std::mt19937 m_gen;
  mutable dist_type m_dist;

  MyFom() : m_gen(rd()), m_dist(0., 1.) { }

  discrete_residual_type createDiscreteTimeResidual() const{
    return discrete_residual_type(N);
  }
  __this_test_phi_type createResultOfDiscreteTimeJacobianActionOn(const __this_test_phi_type & B) const{
    return __this_test_phi_type(N, B.cols());
  }

  template<class StepCountType>
  void discreteTimeResidualAndJacobianAction(StepCountType stepId,
   double time,
   double dt,
   discrete_residual_type & R,
   const __this_test_phi_type & B,
   std::optional<__this_test_phi_type*> /*JA*/,
   const state_type & y_np1,
   const state_type & y_n ) const
  {
    ++count_;
    // std::cout << "FOM: step = " << stepId
    // 	      << " , predictAtTime = " << time
    // 	      << ", dt = " << dt
    // 	      << ", count = " << count_
    // 	      << "\n";

    std::string sp("      ");
    //std::cout << "y_np1, y_n: \n";
    // for (int i=0; i<N; ++i){
    //   std::cout << y_np1(i) << sp << y_n(i) << "\n";
    // }

    if (__this_test_failureHappended){
      const auto romStateStdVec = romStateAtStepStarting.at(stepId);
      auto romState = stdvec_to_eigen(romStateStdVec);
      auto phi = create_gold_phi();
      Eigen::VectorXd tmp = phi * romState;
      EXPECT_TRUE( y_np1.isApprox(y_n) );
      EXPECT_TRUE( y_np1.isApprox(tmp) );
    }
    __this_test_failureHappended = false;

    // if we are one division away from going over the min time step, let is succeed
    // otherwise rnadomly inject failure
    if (dt/reductionFactor > minStepSize){
      const auto coin = m_dist(m_gen);
      //std::cout << coin << "\n";
      if (coin > 0.3){
	throw pressio::eh::DiscreteTimeResidualFailureUnrecoverable();
      }
    }
  }
};

void add_one(Eigen::Matrix<double,-1,1> & a){
  for (int i=0; i<a.size(); ++i){ a(i) += 1.; }
}

struct MyFakeSolver
{
  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    auto R = sys.createResidual();
    auto J = sys.createJacobian();

    for (int k=1; k<=5; ++k){
      try{
	sys.residualAndJacobian(state, R, &J);
      }
      catch (::pressio::eh::ResidualEvaluationFailureUnrecoverable const &e){
	__this_test_failureHappended = true;
	throw ::pressio::eh::NonlinearSolveFailure();
      }
      add_one(state);
      // std::cout << "state = "
      // 		<< state(0) << " "
      // 		<< state(1) << " "
      // 		<< state(2) << "\n";
    }
  }
};
}

TEST(rom_lspg_unsteady, fully_discrete_with_recovery_n2)
{
  /*
    purpose: for the fully discrete lspg API with 2 FOM stencil states
    and time step failure, this test verifies that the FOM states
    passed to the FOM object are always correct, especially in the case
    of failure recovery

    expectation: let's say we are at time step n, then if there is
    a failure/exception during *ANY* iteration of the solve at that step,
    the solve should be attempted again as if nothing happened but with a smaller dt.
    In the case of 2 states, if a failure happens, the FOM states passed
    to the FOM object should be the same becuase the current solution is
    always the initial guess for doing the solve

    some details: for this test, the solver *always* does 5 steps and we fake
    updating the state by simply adding 1 to each element. So not matter what,
    if y_start is the ROM state at the beginning of a time step,
    the solution at the end of the time step should be y_start + 5.
    This is alwasys true. This is just a trick that helps us with this test
    but has no meaning at all. It is just a nmerical trick that helps us verify things.
  */

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

  using fom_t = MyFom;
  using reduced_state_type = Eigen::VectorXd;

  fom_t fomSystem;

  // std::cout << "------------\n";
  // std::cout << "Filling phi \n";
  auto phi = create_gold_phi();
  // std::cout << phi << "\n";
  // std::cout << "------------\n";

  typename fom_t::state_type dummyFomState(N);
  constexpr bool isAffine = false;
  auto space = pressio::rom::create_trial_column_subspace<
    reduced_state_type>(phi, dummyFomState, isAffine);

  auto romState = space.createReducedState();
  romState[0]=0.;
  romState[1]=1.;
  romState[2]=2.;

  auto problem = pressio::rom::lspg::create_unsteady_problem<2>(space, fomSystem);
  auto & stepper = problem.lspgStepper();

  auto dtManager = [](pressio::ode::StepCount /*unused*/,
		      pressio::ode::StepStartAt<double> /*unused*/,
		      pressio::ode::StepSize<double> & dt,
		      pressio::ode::StepSizeMinAllowedValue<double> & minDt,
		      pressio::ode::StepSizeScalingFactor<double> & dtRedFactor)
  {
    dt = 2.;
    minDt = 0.01;
    dtRedFactor=1.2;
  };

  MyFakeSolver solver;
  pressio::ode::advance_to_target_point_with_step_recovery
    (stepper, romState, 0., 2., dtManager, solver);

  PRESSIOLOG_FINALIZE();
}
