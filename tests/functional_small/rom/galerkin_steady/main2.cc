
#include <gtest/gtest.h>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_galerkin_steady.hpp"

namespace{

struct MyFom
{
  using state_type        = Eigen::VectorXd;
  using residual_type     = state_type;
  int N_ = {};

  MyFom(int N): N_(N){}
  residual_type createResidual() const{ return residual_type(N_); }

  Eigen::MatrixXd createResultOfJacobianActionOn(const Eigen::MatrixXd & B) const{
    Eigen::MatrixXd A(N_, B.cols());
    return A;
  }

  void residual(const state_type & u,
		residual_type & r) const
  {}

  template<class OperandType, class ResultType>
  void jacobianAction(const state_type &,
		      OperandType const &,
		      ResultType &) const
  {}

  void residualAndJacobianAction(const state_type & u,
				 residual_type & r,
				 const Eigen::MatrixXd & B,
				 std::optional<Eigen::MatrixXd *> Ain) const
  {
  }

};

}

TEST(rom_galerkin_steady, default_matrix_free)
{
  using namespace pressio;

  /*
    just supposed to compile for now
  */

  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);
  namespace pnonls = pressio::nonlinearsolvers;

  constexpr int N = 8;
  using fom_t = MyFom;
  fom_t fomSystem(N);

  using phi_t = Eigen::Matrix<double, -1,-1>;
  phi_t phi(N, 3);

  using reduced_state_type = Eigen::VectorXd;
  typename fom_t::state_type shift(N);
  auto space = rom::create_trial_column_subspace<reduced_state_type>(phi, shift, false);

  auto problem = rom::galerkin::create_steady_problem(space, fomSystem);

  using tag = linearsolvers::iterative::GMRES;
  auto nonLinSolver = pnonls::experimental::create_matrixfree_newtonkrylov_solver<tag>(problem);

  auto romState = space.createReducedState();
  nonLinSolver.solve(romState);

  PRESSIOLOG_FINALIZE();
}
