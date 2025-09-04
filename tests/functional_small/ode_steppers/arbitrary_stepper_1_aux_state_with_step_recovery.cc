
#include "pressio/ode_steppers.hpp"
#include "pressio/ode_advancers.hpp"

struct MyApp
{
  std::string & checkStr_;
  using independent_variable_type   = double;
  using state_type    = Eigen::VectorXd;
  using discrete_residual_type = state_type;
  using discrete_jacobian_type = Eigen::SparseMatrix<double>;

public:
  MyApp(std::string & strin) : checkStr_(strin){}

  state_type createState() const{ return state_type(3); }

  discrete_residual_type createDiscreteResidual() const{
    discrete_residual_type R(3);
    return R;
  }

  discrete_jacobian_type createDiscreteJacobian() const{
    discrete_jacobian_type J(3,3);
    return J;
  }

  template <typename step_t>
  void discreteResidualAndJacobian(const step_t & step,
				   const independent_variable_type & /*unused*/,
				   const independent_variable_type & dt,
				   discrete_residual_type & R,
				   std::optional<discrete_jacobian_type*> J,
				   const state_type & yn,
				   const state_type & ynm1) const
  {
    std::cout << "yn: " << " "
	      << yn[0] << " "
	      << yn[1] << " "
	      << yn[2] << "\n";

    std::cout << "y_n-1: " << " "
	      << ynm1[0] << " "
	      << ynm1[1] << " "
	      << ynm1[2] << "\n";

    R[0] = ynm1[0]+dt;
    R[1] = ynm1[1]+dt;
    R[2] = ynm1[2]+dt;

    if (step==2){
      if( std::abs(ynm1[0]-1.2) > 1e-13 or
	  std::abs(ynm1[1]-1.4) > 1e-13 or
	  std::abs(ynm1[2]-1.6) > 1e-13){
	checkStr_ = "FAILED";
      }
      if( std::abs(R[0]-1.3) > 1e-13 or
	  std::abs(R[1]-1.5) > 1e-13 or
	  std::abs(R[2]-1.7) > 1e-13){
	checkStr_ = "FAILED";
      }
    }

    if (step==4){
      if( std::abs(ynm1[0]-1.6) > 1e-13 or
    	  std::abs(ynm1[1]-2.2) > 1e-13 or
    	  std::abs(ynm1[2]-2.8) > 1e-13){
    	checkStr_ = "FAILED";
      }
      if( std::abs(R[0]-1.7) > 1e-13 or
    	  std::abs(R[1]-2.3) > 1e-13 or
    	  std::abs(R[2]-2.9) > 1e-13){
    	checkStr_ = "FAILED";
      }
    }

    std::cout << dt << "\n";
    if (step==3 and (dt==0.1 or dt==0.05))
      throw pressio::eh::DiscreteTimeResidualFailureUnrecoverable();

    if (step==3){
      if( std::abs(ynm1[0]-1.4) > 1e-13 or
	  std::abs(ynm1[1]-1.8) > 1e-13 or
	  std::abs(ynm1[2]-2.2) > 1e-13){
	checkStr_ = "FAILED";
      }
      // the only one to succeed at step=3 should be dt=0.025 and such that
      if( std::abs(R[0]-1.425) > 1e-13 or
	  std::abs(R[1]-1.825) > 1e-13 or
	  std::abs(R[2]-2.225) > 1e-13){
	checkStr_ = "FAILED";
      }
    }
  }
};

struct MyFakeSolver
{
  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    // this does not have any meaning, but it mimics the
    // steps happening inside the real pressio nonlin solver

    auto R = sys.createResidual();
    auto J = sys.createJacobian();
    for (int i=0; i<2; ++i)
    {
      std::cout << i << "\n";
      try{
	sys.residualAndJacobian(state, R, &J);
	state(0) += 0.1;
	state(1) += 0.2;
	state(2) += 0.3;
      }
      catch (::pressio::eh::ResidualEvaluationFailureUnrecoverable const &e){
	throw ::pressio::eh::NonlinearSolveFailure();
      }
    }
  }
};


int main()
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::sparse);

  using app_t		= MyApp;
  using state_t	= typename app_t::state_type;

  auto dtManager = [](pressio::ode::StepCount /*unused*/,
		      pressio::ode::StepStartAt<double> /*unused*/,
		      pressio::ode::StepSize<double> & dt,
		      pressio::ode::StepSizeMinAllowedValue<double> & minDt,
		      pressio::ode::StepSizeScalingFactor<double> & dtRedFactor)
  {
    dt = 0.1;
    minDt = 0.01;
    dtRedFactor=2.;
  };

  std::string checkStr= "PASSED";
  app_t appObj(checkStr);
  state_t y(3);
  pressio::ops::fill(y, 1);
  MyFakeSolver solver;

  auto stepperObj = pressio::ode::create_implicit_stepper<2>(appObj);
  pressio::ode::advance_with_step_recovery(stepperObj, y, 0., 0.4, dtManager, solver);

  std::cout << checkStr << std::endl;
  PRESSIOLOG_FINALIZE();
  return 0;
}
