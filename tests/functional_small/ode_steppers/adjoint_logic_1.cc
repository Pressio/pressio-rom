#include <gtest/gtest.h>
#include "pressio/ode_steppers.hpp"
#include "pressio/ode_advancers.hpp"
#include "random"
#include <iomanip>

constexpr int _num_steps = 5;
constexpr int _N = 17;
using _vec_type = Eigen::VectorXd;
using _mat_type = Eigen::MatrixXd;
using _A_collection = std::vector<_mat_type>;

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
  const _A_collection & A_mats_;

public:
  using independent_variable_type = double;
  using CreateBase::state_type;
  using CreateBase::discrete_residual_type;
  using CreateBase::discrete_jacobian_type;

  PrimalApp(const _A_collection & A_mats) : A_mats_(A_mats){}

  void discreteResidualAndJacobian(typename ::pressio::ode::StepCount::value_type stepIn,
				   const independent_variable_type & /*unused*/,
				   const independent_variable_type & dt,
				   discrete_residual_type & R,
				   std::optional<discrete_jacobian_type*> J,
				   const state_type & x_n,
				   const state_type & x_nm1) const
  {
    const int index = stepIn-1;
    std::cout << "primal app: "
	      << " step =  " << stepIn
	      << " index = " << index
	      << std::endl;

    const auto & A = A_mats_[index];
    R = x_n - A * x_nm1;

    if (J){
      J.value()->setIdentity();
    }
  }
};

class AdjointApp : public CreateBase
{
  const _A_collection & A_mats_;

public:
  using independent_variable_type = double;
  using CreateBase::state_type;
  using CreateBase::discrete_residual_type;
  using CreateBase::discrete_jacobian_type;

  AdjointApp(const _A_collection & A_mats) : A_mats_(A_mats){}

  void discreteResidualAndJacobian(typename ::pressio::ode::StepCount::value_type stepIn,
				   const independent_variable_type & /*unused*/,
				   const independent_variable_type & dt,
				   discrete_residual_type & R,
				   std::optional<discrete_jacobian_type*> J,
				   const state_type & w_n,
				   const state_type & w_nm1) const
  {
    const int index = _num_steps - (stepIn-1) - 1;
    std::cout << "adjoint app: "
	      << " step =  " << stepIn
	      << " index = " << index
	      << std::endl;

    const auto & A = A_mats_[index];
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
  }
};

auto create_matrices_for_all_steps(std::mt19937 & gen,
				   std::uniform_real_distribution<double> & dis)
{
  _A_collection As_;
  for (int i=0; i<_num_steps; ++i){
    As_.emplace_back( _mat_type::NullaryExpr(_N,_N, [&](){return dis(gen);}) );
  }
  return As_;
}

TEST(ode, adjoint_test_1)
{
  // see https://github.com/Pressio/pressio/issues/507

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(-1., 1.);

  auto A_mats_ = create_matrices_for_all_steps(gen, dis);

  const auto dtSetter = [](pressio::ode::StepCount /*unused*/,
			   const pressio::ode::StepStartAt<double> & /*unused*/,
			   pressio::ode::StepSize<double> & /*unused*/)
  {
    /* we do not need to set dt because it is not really used in this test */
  };

  //
  // primal
  //
  _vec_type x = _vec_type::NullaryExpr(_N, [&](){return dis(gen);});
  _vec_type x_0 = x;
  std::cout << "x_0 = " << x_0 << std::endl;

  {
    PrimalApp app(A_mats_);
    MyFakeSolver solver;
    auto stepperObj = pressio::ode::create_implicit_stepper<2>(app);
    auto policy = pressio::ode::steps(0., pressio::ode::StepCount(_num_steps), dtSetter);
    pressio::ode::advance(stepperObj, x, policy, solver);
  }

  //
  // adjoint
  //
  _vec_type w = _vec_type::NullaryExpr(_N, [&](){return dis(gen);});
  _vec_type w_0 = w;
  std::cout << "w_0 = " << w_0 << std::endl;

  {
    AdjointApp app(A_mats_);
    MyFakeSolver solver;
    auto stepperObj = pressio::ode::create_implicit_stepper<2>(app);
    auto policy = pressio::ode::steps(0., pressio::ode::StepCount(_num_steps), dtSetter);
    pressio::ode::advance(stepperObj, w, policy, solver);
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
