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
#include "pressio/solvers_linear.hpp"
#include "pressio/solvers_nonlinear.hpp"
#include <Eigen/Core>
#include <optional>

struct MyRosenbrockSystem
{
  using scalar_type   = double;
  using state_type    = Eigen::VectorXd;
  using residual_type = state_type;
  using jacobian_type = Eigen::MatrixXd;

  state_type createState() const{
    auto s = state_type(4);
    s.setZero();
    return s;
  };

  residual_type createResidual() const{
    auto r = residual_type(6);
    r.setZero();
    return r;
  }

  jacobian_type createJacobian() const{
    auto j = jacobian_type(6, 4);
    j.setZero();
    return j;
  }

  void residualAndJacobian(const state_type& x,
			   residual_type& res,
			   std::optional<jacobian_type*> Jin) const
  {
    auto x1 = x(0);
    auto x2 = x(1);
    auto x3 = x(2);
    auto x4 = x(3);

    res(0) = 10.*(x4 - x3*x3);
    res(1) = 10.*(x3 - x2*x2);
    res(2) = 10.*(x2 - x1*x1);
    res(3) = (1.-x1);
    res(4) = (1.-x2);
    res(5) = (1.-x3);

    if (Jin){
      auto & J = *Jin.value();
      J.setZero();
      J(0,2) = -20.*x3; J(0,3) = 10.; // first row
      J(1,1) = -20.*x2; J(1,2) = 10.; // second row
      J(2,0) = -20.*x1; J(2,1) = 10.; // etc
      J(3,0) = -1.;
      J(4,1) = -1.;
      J(5,2) = -1.;
    }
  }
};

int main()
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::sparse);

  namespace pls    = pressio::linearsolvers;
  namespace pnonls = pressio::nonlinearsolvers;

  MyRosenbrockSystem problem;

  using hessian_t = Eigen::MatrixXd;
  using linear_solver_tag = pls::direct::HouseholderQR;
  using linear_solver_t   = pls::Solver<linear_solver_tag, hessian_t>;
  linear_solver_t linearSolver;

  auto x0 = problem.createState();
  x0(0) = 0.0; x0(1) = 1.1; x0(2) = 1.2; x0(3) = 1.5;

  {
    auto solver = pnonls::create_gauss_newton_solver(problem, linearSolver);
    solver.setStopTolerance(1e-5);

    auto x = pressio::ops::clone(x0);
    solver.solve(x);

    // check solution
    std::cout << "Computed solution: \n "
	      << "["
	      << std::setprecision(14)
	      << x(0) << " "
	      << x(1) << " "
	      << x(2) << " "
	      << x(3) << " " << "] \n"
	      << "exact minimum at: \n "
	      << "[1.0, 1.0, 1.0, 1.0] \n";
  }

  {
    auto solver = pnonls::create_levenberg_marquardt_solver(problem, linearSolver);
    solver.setStopTolerance(1e-5);
    solver.setUpdateCriterion(pnonls::Update::LMSchedule2);

    auto x = pressio::ops::clone(x0);
    solver.solve(x);

    // check solution
    std::cout << "Computed solution: \n "
	      << "["
	      << std::setprecision(14)
	      << x(0) << " "
	      << x(1) << " "
	      << x(2) << " "
	      << x(3) << " " << "] \n"
	      << "exact minimum at: \n "
	      << "[1.0, 1.0, 1.0, 1.0] \n";
  }

  std::cout << "PASSED\n";
  PRESSIOLOG_FINALIZE();
  return 0;
}
