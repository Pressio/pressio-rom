/*
//@HEADER
// ************************************************************************
//
// exceptions.hpp
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

#ifndef PRESSIOROM_ODE_EXCEPTIONS_HPP_
#define PRESSIOROM_ODE_EXCEPTIONS_HPP_

#include <exception>

namespace pressio{ namespace eh{

class TimeStepFailure
  : public std::exception
{
  std::string myerr_ = "Time step failed";
  std::string append_ = {};

public:
  TimeStepFailure() = default;
  TimeStepFailure(const TimeStepFailure &) = default;
  TimeStepFailure & operator=(const TimeStepFailure &) = default;
  TimeStepFailure(TimeStepFailure &&) = default;
  TimeStepFailure & operator=(TimeStepFailure &&) = default;
  ~TimeStepFailure() = default;

  explicit TimeStepFailure(std::string append)
    : append_{append}{
    myerr_ += append_;
  }

  virtual const char * what () const noexcept{
    return myerr_.c_str();
   }
};

class VelocityFailureUnrecoverable
  : public std::exception
{
  std::string myerr_ = "Velocity evaluation failed";
  std::string append_ = {};

public:
  VelocityFailureUnrecoverable() = default;
  VelocityFailureUnrecoverable(const VelocityFailureUnrecoverable &) = default;
  VelocityFailureUnrecoverable & operator=(const VelocityFailureUnrecoverable &) = default;
  VelocityFailureUnrecoverable(VelocityFailureUnrecoverable &&) = default;
  VelocityFailureUnrecoverable & operator=(VelocityFailureUnrecoverable &&) = default;
  ~VelocityFailureUnrecoverable() = default;

  explicit VelocityFailureUnrecoverable(std::string append)
    : append_{append}{
    myerr_ += append_;
  }

  virtual const char * what () const noexcept{
    return myerr_.c_str();
   }
};

class DiscreteTimeResidualFailureUnrecoverable
  : public std::exception
{
  std::string myerr_ = "discreteTimeResidual failed";
  std::string append_ = {};

public:
  DiscreteTimeResidualFailureUnrecoverable() = default;
  DiscreteTimeResidualFailureUnrecoverable(const DiscreteTimeResidualFailureUnrecoverable &) = default;
  DiscreteTimeResidualFailureUnrecoverable & operator=(const DiscreteTimeResidualFailureUnrecoverable &) = default;
  DiscreteTimeResidualFailureUnrecoverable(DiscreteTimeResidualFailureUnrecoverable &&) = default;
  DiscreteTimeResidualFailureUnrecoverable & operator=(DiscreteTimeResidualFailureUnrecoverable &&) = default;
  ~DiscreteTimeResidualFailureUnrecoverable() = default;

  explicit DiscreteTimeResidualFailureUnrecoverable(std::string append)
    : append_{append}{
    myerr_ += append_;
  }

  virtual const char * what () const noexcept{
    return myerr_.c_str();
   }
};

}}//end namespace pressio::eh
#endif  // PRESSIOROM_ODE_EXCEPTIONS_HPP_
