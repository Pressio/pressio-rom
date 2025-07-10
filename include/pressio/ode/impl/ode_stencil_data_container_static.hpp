/*
//@HEADER
// ************************************************************************
//
// ode_stencil_data_container_static.hpp
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

#ifndef PRESSIOROM_ODE_IMPL_ODE_STENCIL_DATA_CONTAINER_STATIC_HPP_
#define PRESSIOROM_ODE_IMPL_ODE_STENCIL_DATA_CONTAINER_STATIC_HPP_

namespace pressio{ namespace ode{ namespace impl{

template<class ValueType, std::size_t N, class StencilEndsAtTag>
class StencilDataContainerStaticImpl;

//
// if start at nPlusOne, we have n+1, n, n-1, ...
//
template<class ValueType, std::size_t N>
class StencilDataContainerStaticImpl<ValueType, N, ::pressio::ode::nPlusOne>
{
public:
  using data_type = std::array<ValueType, N>;

private:
  data_type data_;

public:
  template <std::size_t _N = N, std::enable_if_t<_N == 0, int> = 0>
  StencilDataContainerStaticImpl(){}

  // constructor for n == 1
  template <std::size_t _N = N, std::enable_if_t<_N == 1, int> = 0>
  StencilDataContainerStaticImpl(ValueType const & y)
    : data_{::pressio::ops::clone(y)}
  {
    setZero();
  }

  // constructor for n == 2
  template <std::size_t _N = N, std::enable_if_t<_N == 2, int> = 0>
  StencilDataContainerStaticImpl(ValueType const & y)
    : data_{::pressio::ops::clone(y),
            ::pressio::ops::clone(y)}{
    setZero();
  }

  // constructor for n == 3
  template <std::size_t _N = N, std::enable_if_t<_N == 3, int> = 0>
  StencilDataContainerStaticImpl(ValueType const & y)
    : data_{::pressio::ops::clone(y),
            ::pressio::ops::clone(y),
            ::pressio::ops::clone(y)}{
    setZero();
  }

  // constructor for n == 4
  template <std::size_t _N = N, std::enable_if_t<_N == 4, int> = 0>
  StencilDataContainerStaticImpl(ValueType const & y)
    : data_{::pressio::ops::clone(y),
            ::pressio::ops::clone(y),
            ::pressio::ops::clone(y),
            ::pressio::ops::clone(y)}{
    setZero();
  }

  StencilDataContainerStaticImpl(StencilDataContainerStaticImpl const & other) = default;
  StencilDataContainerStaticImpl & operator=(StencilDataContainerStaticImpl const & other) = default;
  StencilDataContainerStaticImpl(StencilDataContainerStaticImpl && other) = default;
  StencilDataContainerStaticImpl & operator=(StencilDataContainerStaticImpl && other) = default;
  ~StencilDataContainerStaticImpl() = default;

public:
  static constexpr std::size_t size(){ return N; }

  ValueType & operator()(::pressio::ode::nPlusOne){
    static_assert( N>=1,
      "Calling operator()(::pressio::ode::nPlusOne) requires N>=1");
    return data_[0];
  }

  ValueType & operator()(::pressio::ode::n){
    static_assert( N>=2,
      "Calling operator()(::pressio::ode::n) requires N>=2");
    return data_[1];
  }

  ValueType & operator()(::pressio::ode::nMinusOne){
    static_assert( N>=3,
      "Calling operator()(::pressio::ode::nMinusOne) requires N>=3");
    return data_[2];
  }

  ValueType & operator()(::pressio::ode::nMinusTwo){
    static_assert( N>=4,
      "Calling operator()(::pressio::ode::nMinusTwo) requires N>=4");
    return data_[3];
  }

  ValueType const & operator()(::pressio::ode::nPlusOne) const {
    static_assert( N>=1,
      "Calling operator()(::pressio::ode::nPlusOne) requires N>=1");
    return data_[0];
  }

  ValueType const & operator()(::pressio::ode::n) const {
    static_assert( N>=2,
      "Calling operator()(::pressio::ode::n) requires N>=2");
    return data_[1];
  }

  ValueType const & operator()(::pressio::ode::nMinusOne) const {
    static_assert( N>=3,
      "Calling operator()(::pressio::ode::nMinusOne) requires N>=3");
    return data_[2];
  }

  ValueType const & operator()(::pressio::ode::nMinusTwo) const {
    static_assert( N>=4,
      "Calling operator()(::pressio::ode::nMinusTwo) requires N>=4");
    return data_[3];
  }

private:
  void setZero(){
    for (auto & it : data_){
      ::pressio::ops::set_zero(it);
    }
  }
};


//
// if we start at n, we have n, n, n-1, ...
//
template<class ValueType, std::size_t N>
class StencilDataContainerStaticImpl<ValueType, N, ::pressio::ode::n>
{
public:
  using data_type = std::array<ValueType, N>;

private:
  data_type data_;

public:
  template <std::size_t _N = N, std::enable_if_t<_N == 0, int> = 0>
  StencilDataContainerStaticImpl(){}

  // constructor for n == 1
  template <std::size_t _N = N, std::enable_if_t<_N == 1, int> = 0>
  StencilDataContainerStaticImpl(ValueType const & y)
    : data_{::pressio::ops::clone(y)}{
    setZero();
  }

  // constructor for n == 2
  template <std::size_t _N = N, std::enable_if_t<_N == 2, int> = 0>
  StencilDataContainerStaticImpl(ValueType const & y)
    : data_{::pressio::ops::clone(y),
            ::pressio::ops::clone(y)}{
    setZero();
  }

  // constructor for n == 3
  template <std::size_t _N = N, std::enable_if_t<_N == 3, int> = 0>
  StencilDataContainerStaticImpl(ValueType const & y)
    : data_{::pressio::ops::clone(y),
            ::pressio::ops::clone(y),
            ::pressio::ops::clone(y)}{
    setZero();
  }

  // constructor for n == 4
  template <std::size_t _N = N, std::enable_if_t<_N == 4, int> = 0>
  StencilDataContainerStaticImpl(ValueType const & y)
    : data_{::pressio::ops::clone(y),
            ::pressio::ops::clone(y),
            ::pressio::ops::clone(y),
            ::pressio::ops::clone(y)}{
    setZero();
  }

  StencilDataContainerStaticImpl(StencilDataContainerStaticImpl const & other) = default;
  StencilDataContainerStaticImpl & operator=(StencilDataContainerStaticImpl const & other) = default;
  StencilDataContainerStaticImpl(StencilDataContainerStaticImpl && other) = default;
  StencilDataContainerStaticImpl & operator=(StencilDataContainerStaticImpl && other) = default;
  ~StencilDataContainerStaticImpl() = default;

public:
  static constexpr std::size_t size(){ return N; }

  ValueType & operator()(::pressio::ode::n){
    static_assert( N>=1,
      "Calling operator()(::pressio::ode::n) requires N>=1");
    return data_[0];
  }

  ValueType & operator()(::pressio::ode::nMinusOne){
    static_assert( N>=2,
      "Calling operator()(::pressio::ode::nMinusOne) requires N>=2");
    return data_[1];
  }

  ValueType & operator()(::pressio::ode::nMinusTwo){
    static_assert( N>=3,
      "Calling operator()(::pressio::ode::nMinusTwo) requires N>=3");
    return data_[2];
  }

  ValueType & operator()(::pressio::ode::nMinusThree){
    static_assert( N>=4,
      "Calling operator()(::pressio::ode::nMinusThree) requires N>=4");
    return data_[3];
  }

  ValueType const & operator()(::pressio::ode::n) const {
    static_assert( N>=1,
      "Calling operator()(::pressio::ode::n) requires N>=1");
    return data_[0];
  }

  ValueType const & operator()(::pressio::ode::nMinusOne) const {
    static_assert( N>=2,
      "Calling operator()(::pressio::ode::nMinusOne) requires N>=2");
    return data_[1];
  }

  ValueType const & operator()(::pressio::ode::nMinusTwo) const {
    static_assert( N>=3,
      "Calling operator()(::pressio::ode::nMinusTwo) requires N>=3");
    return data_[2];
  }

  ValueType const & operator()(::pressio::ode::nMinusThree) const {
    static_assert( N>=4,
      "Calling operator()(::pressio::ode::nMinusThree) requires N>=4");
    return data_[3];
  }

private:
  void setZero(){
    for (auto & it : data_){
      ::pressio::ops::set_zero(it);
    }
  }
};

}}}
#endif  // PRESSIOROM_ODE_IMPL_ODE_STENCIL_DATA_CONTAINER_STATIC_HPP_
