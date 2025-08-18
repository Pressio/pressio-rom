#include <gtest/gtest.h>
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "random"
#include <iomanip>

constexpr int _num_steps = 5;
constexpr int _N = 17;
using _vec_type = Eigen::VectorXd;
using _mat_type = Eigen::MatrixXd;
using _B_collection = std::vector<_mat_type>;

class CreateBase{
public:
  using state_type = _vec_type;
  using discrete_residual_type = _vec_type;
  using discrete_jacobian_type = _mat_type;

  state_type createState() const{
    state_type a(_N);
    a.setZero();
    return a;
  }

  discrete_residual_type createDiscreteResidual() const{
    discrete_residual_type a(_N);
    a.setZero();
    return a;
  }

  discrete_jacobian_type createDiscreteJacobian() const{
    discrete_jacobian_type J(_N, _N);
    return J;
  }
};

class PrimalApp : public CreateBase
{
  const _B_collection & B_mats_;
  const _vec_type & stepSizes_;

public:
  using independent_variable_type = double;
  using CreateBase::state_type;
  using CreateBase::discrete_residual_type;
  using CreateBase::discrete_jacobian_type;

  PrimalApp(const _B_collection & mats,
	    const _vec_type & stepSizes)
    : B_mats_(mats), stepSizes_(stepSizes){}

  void discreteResidualAndJacobian(typename ::pressio::ode::StepCount::value_type stepIn,
				   const independent_variable_type & /*unused*/,
				   const independent_variable_type & dt,
				   discrete_residual_type & R,
				   std::optional<discrete_jacobian_type*> J,
				   const state_type & x_n,
				   const state_type & x_nm1) const
  {
    const int index = stepIn-1;
    EXPECT_DOUBLE_EQ(dt, stepSizes_[index]);

    std::cout << "primal app: "
	      << " step =  " << stepIn
	      << " index = " << index
	      << " dt = " << dt
	      << std::endl;

    auto A = B_mats_[index];
    ::pressio::ops::add_to_diagonal(A, dt);
    R = x_n - A * x_nm1;
    if (J){
      J.value()->setIdentity();
    }
  }
};

class AdjointApp : public CreateBase
{
  const _B_collection & B_mats_;
  const _vec_type & stepSizes_;

public:
  using independent_variable_type = double;
  using CreateBase::state_type;
  using CreateBase::discrete_residual_type;
  using CreateBase::discrete_jacobian_type;

  AdjointApp(const _B_collection & mats,
	     const _vec_type & stepSizes)
    : B_mats_(mats), stepSizes_(stepSizes){}

  void discreteResidualAndJacobian(typename ::pressio::ode::StepCount::value_type stepIn,
				   const independent_variable_type & /*unused*/,
				   const independent_variable_type & dt,
				   discrete_residual_type & R,
				   std::optional<discrete_jacobian_type*> J,
				   const state_type & w_n,
				   const state_type & w_nm1) const
  {
    const int index = _num_steps - (stepIn-1) - 1;

    EXPECT_TRUE(dt < 0.);
    // use negative because adjoint goes backward
    const auto dtToUse = -dt;

    EXPECT_DOUBLE_EQ(dt, stepSizes_[index]);

    std::cout << "adjoint app: "
	      << " step =  " << stepIn
	      << " index = " << index
	      << " dt = " << dt
	      << std::endl;

    auto A = B_mats_[index];
    ::pressio::ops::add_to_diagonal(A, dtToUse);
    R = w_n - A.transpose() * w_nm1;
    if (J){
      J.value()->setIdentity();
    }
  }
};

struct MyFakeSolver{
  template<class system_t, class state_t>
  void solve(const system_t & sys, state_t & state)
  {
    auto R = sys.createResidual();
    auto J = sys.createJacobian();

    sys.residualAndJacobian(state, R, &J);
    assert(J == _mat_type::Identity(J.rows(), J.cols()));
    Eigen::VectorXd correction = R;
    state -= correction;
    //std::cout << "solve: state = " << state << std::endl;
  }
};

auto create_matrices_for_all_steps(std::mt19937 & gen,
				   std::uniform_real_distribution<double> & dis)
{
  _B_collection As_;
  for (int i=0; i<_num_steps; ++i){
    As_.emplace_back( _mat_type::NullaryExpr(_N,_N, [&](){return dis(gen);}) );
  }
  return As_;
}

struct PrimalObserver
{
  double finalTime_;
  _vec_type stepSizes_;

  PrimalObserver(double finalTime, const _vec_type & stepSizes)
    : finalTime_(finalTime), stepSizes_(stepSizes){}

  void operator()(pressio::ode::StepCount stepIn,
		  double timeIn,
		  const Eigen::VectorXd & /**/) const
  {
    std::cout << stepIn.get() << " " << timeIn << std::endl;
    if (stepIn.get() == 0){ EXPECT_DOUBLE_EQ( timeIn, 0.); }
    else { EXPECT_DOUBLE_EQ( timeIn, stepSizes_.head(stepIn.get()).sum()); }
  }
};

struct AdjointObserver
{
  double finalTime_;
  _vec_type stepSizes_;
  AdjointObserver(double finalTime,
		  const _vec_type & stepSizes)
    : finalTime_(finalTime), stepSizes_(stepSizes){}

  void operator()(pressio::ode::StepCount stepIn,
		  double timeIn,
		  const Eigen::VectorXd & /**/) const
  {
    std::cout << stepIn.get() << " " << timeIn << std::endl;
    if (stepIn.get() == 0){
      EXPECT_DOUBLE_EQ( timeIn, finalTime_); }
    else {
      EXPECT_DOUBLE_EQ( timeIn, finalTime_ - std::abs(stepSizes_.tail(stepIn.get()).sum()));
    }
  }
};

TEST(ode, adjoint_test_2)
{
  // see https://github.com/Pressio/pressio/issues/507

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(-1., 1.);
  auto B_mats_ = create_matrices_for_all_steps(gen, dis);

  _vec_type stepSizes(_num_steps);
  stepSizes << 0.5, 0.3, 1.2, 0.7, 0.25;
  const auto finalTime = stepSizes.sum();
  assert( stepSizes.size() == _num_steps);

  //
  // primal
  //
  _vec_type x = _vec_type::NullaryExpr(_N, [&](){return dis(gen);});
  _vec_type x_0 = x;
  std::cout << "x_0 = " << x_0 << std::endl;

  {
    const auto dtSetter = [=](pressio::ode::StepCount stepIn,
			      const pressio::ode::StepStartAt<double> & /*unused*/,
			      pressio::ode::StepSize<double> & dt)
    {
      // need to use step - 1 because pressio steps are
      // enumerated as 1,2,3...
      dt = stepSizes[stepIn.get()-1];
    };

    PrimalApp app(B_mats_, stepSizes);
    MyFakeSolver solver;
    PrimalObserver obs(finalTime, stepSizes);
    auto stepperObj = pressio::ode::create_implicit_stepper<2>(app);
    auto policy = pressio::ode::steps(0., pressio::ode::StepCount(_num_steps), dtSetter);
    pressio::ode::advance(stepperObj, x, policy, solver, obs);
  }

  //
  // adjoint
  //
  _vec_type w = _vec_type::NullaryExpr(_N, [&](){return dis(gen);});
  _vec_type w_0 = w;
  std::cout << "w_0 = " << w_0 << std::endl;

  {
    pressio::ops::scale(stepSizes, -1.0);
    const auto dtSetter = [=](pressio::ode::StepCount stepIn,
			      const pressio::ode::StepStartAt<double> & /*unused*/,
			      pressio::ode::StepSize<double> & dt)
    {
      // need to use step - 1 because pressio steps are
      // enumerated as 1,2,3...
      const int index = _num_steps - (stepIn.get()-1) - 1;
      dt = stepSizes[index];
    };

    AdjointApp app(B_mats_, stepSizes);
    MyFakeSolver solver;
    AdjointObserver obs(finalTime, stepSizes);
    auto stepperObj = pressio::ode::create_implicit_stepper<2>(app);
    auto policy = pressio::ode::steps(finalTime, pressio::ode::StepCount(_num_steps), dtSetter);
    pressio::ode::advance(stepperObj, w, policy, solver, obs);
  }

  const auto lhs = pressio::ops::dot(w, x_0);
  const auto rhs = pressio::ops::dot(x, w_0);
  std::cout << " dot(w_n , x_0) = " << std::setprecision(14) << lhs << std::endl;
  std::cout << " dot(x_n , w_0) = " << std::setprecision(14) << rhs << std::endl;
  std::cout << lhs - rhs << std::endl;

  // this tolerance cannot be too small because all A matrices contain
  // elements that are not small so when doing all products we get large
  // large numbers at the end
  ASSERT_NEAR(lhs, rhs, 1e-12);
}
