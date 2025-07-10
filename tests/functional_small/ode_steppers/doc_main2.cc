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
#include "pressio/ode_steppers_explicit.hpp"
#include <Eigen/Core>

template<class ScalarType>
class Lorenz3{
  const int N_ = 3;
  const ScalarType rho_   = static_cast<ScalarType>(28);
  const ScalarType sigma_ = static_cast<ScalarType>(10);
  const ScalarType beta_  = static_cast<ScalarType>(8)/3;

public:
  using independent_variable_type = ScalarType;
  using state_type                = Eigen::Matrix<ScalarType,-1,1>;
  using rhs_type      = state_type;

  state_type createState() const{
    auto s = state_type(N_);
    s.setConstant(0);
    return s;
  };

  rhs_type createRhs() const{
    auto v = rhs_type(N_);
    v.setConstant(0);
    return v;
  };

  void rhs(const state_type & state,
	   const independent_variable_type timeIn,
	   rhs_type & rhs) const
  {
    const auto x = state(0);
    const auto y = state(1);
    const auto z = state(2);

    rhs(0) = sigma_ * (y - x);
    rhs(1) = x * (rho_ - z) - y;
    rhs(2) = x * y - beta_ * z;
  };
};

template<class ScalarType>
class StateObserver
{
public:
  explicit StateObserver(int sampleFreq)
    : myfile_("state_snapshots.bin",  std::ios::out | std::ios::binary),
      sampleFreq_(sampleFreq){}

  ~StateObserver(){ myfile_.close(); }

  template<typename TimeType>
  void operator()(pressio::ode::StepCount step,
		  const TimeType /*timeIn*/,
		  const Eigen::Matrix<ScalarType,-1,1> & state)
  {
    if (step.get() % sampleFreq_ == 0){
      const std::size_t ext = state.size()*sizeof(ScalarType);
      myfile_.write(reinterpret_cast<const char*>(&state(0)), ext);
    }
  }

private:
  std::ofstream myfile_;
  const int sampleFreq_ = {};
};

int main(int argc, char *argv[])
{
  PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::sparse);

  using scalar_type = double;
  Lorenz3<scalar_type> problem;

  namespace pode = pressio::ode;
  constexpr auto scheme = pode::StepScheme::RungeKutta4;
  auto stepperObj = pode::create_explicit_stepper(scheme, problem);

  auto y = problem.createState();
  y(0) = 0.5; y(1) = 1.0; y(2) = 1.;
  const scalar_type startTime{0.0};
  const scalar_type finalTime{40.0};
  const int stateSampleFreq = 2;
  pode::advance_to_target_point(stepperObj, y, startTime,
				finalTime,
				  [=](const pode::StepCount & /*unused*/,
				     const pode::StepStartAt<scalar_type> & /*unused*/,
				     pode::StepSize<scalar_type> & dt)
				  { dt = 0.01; },
				StateObserver<scalar_type>{stateSampleFreq});

  std::cout << "PASSED\n";
  PRESSIOLOG_FINALIZE();
  return 0;
}
