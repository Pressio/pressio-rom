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
#include "pressio/ode_advancers.hpp"
#include "pressio/ode_steppers_explicit.hpp"
#include <Eigen/Core>

template<class ScalarType>
class SimpleSystem{
public:
  using independent_variable_type = ScalarType;
  using state_type                = Eigen::Matrix<ScalarType,-1,1>;
  using rhs_type	          = state_type;

  state_type createState() const{
    auto s = state_type(3);
    s.setConstant(0);
    return s;
  };

  rhs_type createRhs() const{
    auto v = rhs_type(3);
    v.setConstant(0);
    return v;
  };

  void rhs(const state_type & state,
	   const independent_variable_type timeIn,
	   rhs_type & rhs) const
  {
    constexpr ScalarType ten{10};
    rhs(0) = ten * state(0);
    rhs(1) = ten * state(1);
    rhs(2) = ten * state(2);
  };
};

int main(int argc, char *argv[])
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::sparse);

  using scalar_type = double;
  SimpleSystem<scalar_type> problem;

  auto y = problem.createState();

  namespace pode = pressio::ode;
  constexpr auto scheme = pode::StepScheme::ForwardEuler;
  auto stepperObj = pode::create_explicit_stepper(scheme, problem);

  {
    /*
      run a fixed number of steps
    */
    y(0) = 1.; y(1) = 2.; y(2) = 3.;
    const scalar_type dt{0.1};
    const scalar_type startTime{0.0};
    pode::advance_n_steps(stepperObj, y, startTime, dt, pode::StepCount{1});
    // check solution
    std::cout << "Computed solution: ["
	      << y(0) << " " << y(1) << " " << y(2) << "] "
	      << "Expected solution: [2,4,6] \n";
  }
  std::cout << std::endl;

  {
    /*
      run until target time
    */
    y(0) = 1.; y(1) = 2.; y(2) = 3.;
    const scalar_type startTime{0.0};
    const scalar_type finalTime{1.0};
    pode::advance_to_target_point(stepperObj, y, startTime, finalTime,
				  [](const pode::StepCount & /*unused*/,
				     const pode::StepStartAt<scalar_type> & /*unused*/,
				     pode::StepSize<scalar_type> & dt)
				  {
				    dt = 0.1;
				  });
    // check solution
    std::cout << "Computed solution: ["
	      << y(0) << " " << y(1) << " " << y(2) << "] "
	      << "Expected solution: [1024, 2048, 3074] \n";
  }

  std::cout << "PASSED\n";
  PRESSIOLOG_FINALIZE();
  return 0;
}
