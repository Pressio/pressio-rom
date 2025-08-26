/*
//@HEADER
// ************************************************************************
//
// main1.cc
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
#include "pressio/solvers_nonlinear.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <optional>

struct MyProblem
{
  using state_type    = Eigen::VectorXd;
  using residual_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

  state_type createState() const {
    state_type a(2);
    a.setZero();
    return a;
  }

  residual_type createResidual() const {
    residual_type a(2);
    a.setZero();
    return a;
  }

  jacobian_type createJacobian() const {
    jacobian_type a(2, 2);
    a.setZero();
    return a;
  }

  void residualAndJacobian(const state_type& x,
			   residual_type& res,
			   std::optional<jacobian_type*> Jin) const
  {
    res(0) = x(0)*x(0) + x(1)*x(1) - 4.0;
    res(1) = x(0)*x(1) - 1.0;

    if (Jin){
     auto & jac = *Jin.value();
      jac.coeffRef(0, 0) = 2*x(0);
      jac.coeffRef(0, 1) = 2*x(1);
      // Have incorrect entries so that line search is required
      jac.coeffRef(1, 0) = 0.1*x(1);
      jac.coeffRef(1, 1) = 0.1*x(0);
    }
  }
};

class MyLinearSolver{
  using A_t = typename MyProblem::jacobian_type;
  using b_t = typename MyProblem::residual_type;
  using x_t = typename MyProblem::state_type;

public:
  void solve(A_t const & A, b_t const & b, x_t & x)
  {
    Eigen::BiCGSTAB< A_t > solver;
    solver.compute(A);
    x = solver.solve(b);
  }
};

int main()
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::sparse);

  namespace pnonls = pressio::nlsol;

  MyProblem problem;
  MyLinearSolver linSolver;
  auto solver = pnonls::create_newton_solver(problem, linSolver);
  const auto updateMethod = pnonls::Update::BacktrackStrictlyDecreasingObjective;
  solver.setUpdateCriterion(updateMethod);

  auto x = problem.createState();
  x(0) = 0.3; x(1) = 0.4;
  solver.solve(x);

  // check solution
  std::cout << "Computed solution: \n "
	    << "["
	    << std::setprecision(14)
	    << x(0) << " "
	    << x(1) << " " << "] \n"
	    << "Expected solution: \n "
	    << "[1.93185165, 0.51763809] \n";

  std::cout << "PASSED\n";
  PRESSIOLOG_FINALIZE();
  return 0;
}
