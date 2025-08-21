
#include <gtest/gtest.h>
#include "pressio/ode_concepts.hpp"

struct FakeIndVarTypeForTesting{
  operator double(){ return double{}; }
};

//
// rhs only
//
struct System1{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = Eigen::VectorXd;
  using rhs_type = Eigen::VectorXd;

  state_type createState() const{ return state_type(); }
  rhs_type createRhs() const{ return rhs_type(); }

  void rhs(const state_type &        /*unused*/,
	   independent_variable_type /*unused*/,
	   rhs_type &    /*unused*/) const{}
};

struct System2{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = Eigen::VectorXd;
  using rhs_type = Eigen::VectorXd;

  // omit createRhs so that it is invalid

  void rhs(const state_type &        /*unused*/,
	   independent_variable_type /*unused*/,
	   rhs_type &    /*unused*/) const{}
};

struct System3{
  // omit independent_variable_type so that it is valid
  using state_type = Eigen::VectorXd;
  using rhs_type = Eigen::VectorXd;
  rhs_type createRhs() const{ return rhs_type(); }

  void rhs(const state_type &       /*unused*/,
	   FakeIndVarTypeForTesting /*unused*/,
	   rhs_type &   /*unused*/) const{}
};

//
// rhs and mass matrix
//
struct System4{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = Eigen::VectorXd;
  using rhs_type = Eigen::VectorXd;
  using mass_matrix_type = Eigen::SparseMatrix<double>;

  state_type createState() const{ return state_type(); }
  rhs_type createRhs() const{ return rhs_type(); }
  mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }

  void massMatrixAndRhs(const state_type &        /*unused*/,
			independent_variable_type /*unused*/,
			mass_matrix_type &        /*unused*/,
			rhs_type &    /*unused*/) const{}
};

//
// rhs and jacobian
//
struct System6{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = Eigen::VectorXd;
  using rhs_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;

  state_type createState() const{ return state_type(); }
  rhs_type createRhs() const{ return rhs_type(); }
  jacobian_type createJacobian() const{ return jacobian_type(); }

  void rhsAndJacobian(const state_type &        /*unused*/,
		      independent_variable_type /*unused*/,
		      rhs_type &    /*unused*/,
		      std::optional<jacobian_type*> /*unused*/) const{}
};

//
// rhs and jacobian and mass matrix
//
struct System7{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = Eigen::VectorXd;
  using rhs_type = Eigen::VectorXd;
  using jacobian_type = Eigen::MatrixXd;
  using mass_matrix_type = Eigen::SparseMatrix<double>;

  state_type createState() const{ return state_type(); }
  rhs_type createRhs() const{ return rhs_type(); }
  jacobian_type createJacobian() const{ return jacobian_type(); }
  mass_matrix_type createMassMatrix() const{ return mass_matrix_type(); }

  void massMatrixAndRhsAndJacobian(const state_type &        /*unused*/,
				   independent_variable_type /*unused*/,
				   mass_matrix_type &        /*unused*/,
				   rhs_type &    /*unused*/,
				   std::optional<jacobian_type*> /*unused*/) const{}
};

//
// fully discrete
//
struct System9{
  using independent_variable_type = FakeIndVarTypeForTesting;
  using state_type = Eigen::VectorXd;
  using discrete_residual_type = Eigen::VectorXd;
  using discrete_jacobian_type = Eigen::MatrixXd;

  state_type createState() const{ return state_type(); }

  discrete_residual_type createDiscreteResidual() const{
    return discrete_residual_type(); }
  discrete_jacobian_type createDiscreteJacobian() const{
    return discrete_jacobian_type(); }

  template<class StepCountType>
  void discreteResidualAndJacobian(StepCountType,
				   independent_variable_type /*unused*/,
				   independent_variable_type /*unused*/,
				   discrete_residual_type &  /*unused*/,
           std::optional<discrete_jacobian_type*> /*unused*/,
				   const state_type &        /*unused*/) const{}

  template<class StepCountType>
  void discreteResidualAndJacobian(StepCountType,
				   independent_variable_type /*unused*/,
				   independent_variable_type /*unused*/,
				   discrete_residual_type &  /*unused*/,
           std::optional<discrete_jacobian_type*> /*unused*/,
				   const state_type &        /*unused*/,
				   const state_type &        /*unused*/) const{}
};

#include <gtest/gtest.h>

TEST(ode, concepts)
{
  using namespace pressio::ode;

  static_assert( PRESSIO_VALUE_OF(OdeSystem<System1>), "");
  static_assert(!PRESSIO_VALUE_OF(CompleteOdeSystem<System1>), "");

  static_assert(!PRESSIO_VALUE_OF(OdeSystem<System2>), "");
  static_assert(!PRESSIO_VALUE_OF(CompleteOdeSystem<System2>), "");

  static_assert(!PRESSIO_VALUE_OF(OdeSystem<System3>), "");
  static_assert(!PRESSIO_VALUE_OF(CompleteOdeSystem<System3>), "");

  static_assert( PRESSIO_VALUE_OF(OdeSystemFusingMassMatrixAndRhs<System4>), "");

  static_assert(!PRESSIO_VALUE_OF(CompleteOdeSystem<System6>), "");
  static_assert( PRESSIO_VALUE_OF(CompleteOdeSystem<System7>), "");

  static_assert( PRESSIO_VALUE_OF((FullyDiscreteSystemWithJacobian<System9, 1>)), "");
  static_assert( PRESSIO_VALUE_OF((FullyDiscreteSystemWithJacobian<System9, 2>)), "");
  static_assert(!PRESSIO_VALUE_OF((FullyDiscreteSystemWithJacobian<System9, 3>)), "");
  static_assert(!PRESSIO_VALUE_OF(OdeSystem<System9>), "");
  static_assert(!PRESSIO_VALUE_OF(CompleteOdeSystem<System9>), "");
}
