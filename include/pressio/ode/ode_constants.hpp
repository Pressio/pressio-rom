/*
//@HEADER
// ************************************************************************
//
// ode_constants.hpp
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

#ifndef PRESSIOROM_ODE_ODE_CONSTANTS_HPP_
#define PRESSIOROM_ODE_ODE_CONSTANTS_HPP_

namespace pressio{ namespace ode{

constexpr typename StepCount::value_type first_step_value = 1;

namespace constants{

template <typename scalar_t>
struct Constants
{
  static constexpr scalar_t negOne(){ return static_cast<scalar_t>(-1); }
  static constexpr scalar_t zero()  { return static_cast<scalar_t>(0);  }
  static constexpr scalar_t one()   { return static_cast<scalar_t>(1);  }
  static constexpr scalar_t two()   { return static_cast<scalar_t>(2);  }
  static constexpr scalar_t three() { return static_cast<scalar_t>(3);  }
  static constexpr scalar_t four()  { return static_cast<scalar_t>(4);  }
  static constexpr scalar_t six()   { return static_cast<scalar_t>(6);  }

  static constexpr scalar_t negOneHalf() { return negOne()/two(); }
  static constexpr scalar_t oneOvThree()  { return one()/three(); }
  static constexpr scalar_t twoOvThree()  { return two()/three(); }
  static constexpr scalar_t threeOvFour() { return three()/four(); }
  static constexpr scalar_t fourOvThree() { return four()/three(); }
  static constexpr scalar_t threeOvTwo()  { return three()/two(); }
  static constexpr scalar_t fourInv()     { return one()/four(); }
};

template <typename scalar_t>
struct bdf1{
  using cnst = Constants<scalar_t>;
  static constexpr scalar_t c_np1_= cnst::one();
  static constexpr scalar_t c_n_  = cnst::negOne();
  static constexpr scalar_t c_f_  = cnst::negOne();
};

template <typename scalar_t>
struct bdf2{
  using cnst = Constants<scalar_t>;
  static constexpr scalar_t c_np1_ = cnst::one();
  static constexpr scalar_t c_n_   = cnst::negOne()*cnst::fourOvThree();
  static constexpr scalar_t c_nm1_ = cnst::oneOvThree();
  static constexpr scalar_t c_f_   = cnst::negOne()*cnst::twoOvThree();
};

template <typename scalar_t>
struct cranknicolson{
  using cnst = Constants<scalar_t>;
  static constexpr scalar_t c_np1_  = cnst::one();
  static constexpr scalar_t c_n_    = cnst::negOne();
  static constexpr scalar_t c_fnp1_ = cnst::negOneHalf();
  static constexpr scalar_t c_fn_   = cnst::negOneHalf();
};

}}}//end namespace pressio::ode::constants

#endif  // PRESSIOROM_ODE_ODE_CONSTANTS_HPP_
