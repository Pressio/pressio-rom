/*
//@HEADER
// ************************************************************************
//
// main.cc
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#include "pressio/ode_advancers.hpp"
#include "pressio/ode_steppers.hpp"
#include <Eigen/Core>

template<class ScalarType>
class Lorenz3{
  const int N_ = 3;
  const ScalarType rho_   = static_cast<ScalarType>(28);
  const ScalarType sigma_ = static_cast<ScalarType>(10);
  const ScalarType beta_  = static_cast<ScalarType>(8)/3;

  using trip_t = Eigen::Triplet<ScalarType>;
  mutable std::vector<trip_t> tripletList_;

public:
  using independent_variable_type = ScalarType;
  using state_type                = Eigen::Matrix<ScalarType,-1,1>;
  using rhs_type		  = state_type;
  // since the Jacobain of Lorenz3 is fairly dens,
  // we could use a dense matrix for it, but we use a sparse one
  // here for demonstration purposes
  using jacobian_type             = Eigen::SparseMatrix<ScalarType>;

  Lorenz3(){
    allocateTriplets();
  }

  state_type createState() const{
    state_type s(N_);
    s.setConstant(0);
    return s;
  };

  rhs_type createRhs() const{
    rhs_type v(N_);
    v.setConstant(0);
    return v;
  };

  jacobian_type createJacobian() const{
    jacobian_type J(N_, N_);
    J.setFromTriplets(tripletList_.begin(), tripletList_.end());
    return J;
  };

  void rhsAndJacobian(const state_type & state,
		      const independent_variable_type timeIn,
		      rhs_type & rhs,
		      std::optional<jacobian_type*> Jopt) const
  {
    const auto x = state(0);
    const auto y = state(1);
    const auto z = state(2);

    rhs(0) = sigma_ * (y - x);
    rhs(1) = x * (rho_ - z) - y;
    rhs(2) = x * y - beta_ * z;

    if (Jopt)
    {
      auto & J = *Jopt.value();
      tripletList_[0] = {0, 0, -sigma_};
      tripletList_[1] = {0, 1,  sigma_};

      tripletList_[2] = {1, 0, rho_-z};
      tripletList_[3] = {1, 1, -1.};
      tripletList_[4] = {1, 2, -x};

      tripletList_[5] = {2, 0, y};
      tripletList_[6] = {2, 1, x};
      tripletList_[7] = {2, 2, -beta_};
      J.setFromTriplets(tripletList_.begin(), tripletList_.end());
    }
  };

private:
  void allocateTriplets()
  {
    // first row
    tripletList_.push_back( trip_t( 0, 0, 0.) );
    tripletList_.push_back( trip_t( 0, 1, 0.) );
    // second and third row
    for (int i=1; i<=2; ++i){
      tripletList_.push_back( trip_t( i, 0, 0.) );
      tripletList_.push_back( trip_t( i, 1, 0.) );
      tripletList_.push_back( trip_t( i, 2, 0.) );
    }
  }
};

int main(int argc, char *argv[])
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::sparse);

  using scalar_type = double;
  using problem_type = Lorenz3<scalar_type>;
  problem_type problem;

  namespace pode = pressio::ode;
  namespace plins  = pressio::linearsolvers;
  namespace pnonls = pressio::nonlinearsolvers;

  // create stepper
  constexpr auto scheme = pode::StepScheme::BDF2;
  auto stepper = pode::create_implicit_stepper(scheme, problem);

  // linear and nonlinear solver
  using fom_jacobian_t = typename problem_type::jacobian_type;
  using linear_solver_t = plins::Solver<plins::iterative::Bicgstab, fom_jacobian_t>;
  linear_solver_t linearSolver;

  auto nonLinearSolver = pnonls::create_newton_solver(stepper, linearSolver);
  nonLinearSolver.setStopTolerance(1e-13);

  auto y = problem.createState();
  y(0) = 0.5; y(1) = 1.0; y(2) = 1.;
  const scalar_type startTime{0.0};
  const scalar_type finalTime{40.0};
  auto policy = pode::to_time(startTime, finalTime,
			      [=](const pode::StepCount & /*unused*/,
				  const pode::StepStartAt<scalar_type> & /*unused*/,
				  pode::StepSize<scalar_type> & dt)
			      { dt = 0.01; });

  std::ofstream myfile("state_snapshots.bin",  std::ios::out | std::ios::binary);
  pode::advance(stepper, y, policy, nonLinearSolver,
		// lambda to observe and write to file the state
		[&](pressio::ode::StepCount step,
		    const scalar_type /*timeIn*/,
		    const typename problem_type::state_type & state)
		{
		  std::cout << "saving state at = " << step.get() << "\n";
		  const std::size_t ext = state.size()*sizeof(scalar_type);
		  myfile.write(reinterpret_cast<const char*>(&state(0)), ext);
		});
  myfile.close();

  std::cout << "PASSED\n";
  PRESSIOLOG_FINALIZE();
  return 0;
}
