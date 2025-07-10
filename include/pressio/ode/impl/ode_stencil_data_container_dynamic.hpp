/*
//@HEADER
// ************************************************************************
//
// ode_stencil_data_container_dynamic.hpp
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

#ifndef PRESSIOROM_ODE_IMPL_ODE_STENCIL_DATA_CONTAINER_DYNAMIC_HPP_
#define PRESSIOROM_ODE_IMPL_ODE_STENCIL_DATA_CONTAINER_DYNAMIC_HPP_

namespace pressio{ namespace ode{ namespace impl{

template<class ValueType, class StencilEndsAtTag>
class StencilDataContainerDynImpl;

//
// specialize for nPlusOne:
// so we have n+1, n, n-1, ...
//
template<class ValueType>
class StencilDataContainerDynImpl<ValueType, ::pressio::ode::nPlusOne>
{
public:
  using data_type = std::vector<ValueType>;

public:
  StencilDataContainerDynImpl() = default;

  StencilDataContainerDynImpl(std::initializer_list<ValueType> il)
    : data_{il}, size_(data_.size())
  {
    for (auto & it : data_){
      ::pressio::ops::set_zero(it);
    }
  }

  StencilDataContainerDynImpl(StencilDataContainerDynImpl const & other) = default;
  StencilDataContainerDynImpl & operator=(StencilDataContainerDynImpl const & other) = default;
  StencilDataContainerDynImpl(StencilDataContainerDynImpl && other) = default;
  StencilDataContainerDynImpl & operator=(StencilDataContainerDynImpl && other) = default;
  ~StencilDataContainerDynImpl() = default;

public:
  std::size_t size() const{ return size_; }

  // non-const overloads
  ValueType & operator()(::pressio::ode::nPlusOne){
    assert(size_>=1);
    return data_[0];
  }

  ValueType & operator()(::pressio::ode::n){
    assert( size_>=2);
    return data_[1];
  }

  ValueType & operator()(::pressio::ode::nMinusOne){
    assert( size_>=3);
    return data_[2];
  }

  ValueType & operator()(::pressio::ode::nMinusTwo){
    assert( size_>=4);
    return data_[3];
  }

  // const overloads
  ValueType const & operator()(::pressio::ode::nPlusOne) const {
    assert( size_>=1);
    return data_[0];
  }

  ValueType const & operator()(::pressio::ode::n) const {
    assert( size_>=2);
    return data_[1];
  }

  ValueType const & operator()(::pressio::ode::nMinusOne) const {
    assert( size_>=3);
    return data_[2];
  }

  ValueType const & operator()(::pressio::ode::nMinusTwo) const {
    assert( size_>=4);
    return data_[3];
  }

private:
  data_type data_;
  std::size_t size_;
};


//
// specialize for n:
// so we have n, n-1, n-2, ...
//
template<class ValueType>
class StencilDataContainerDynImpl<ValueType, ::pressio::ode::n>
{
public:
  using data_type = std::vector<ValueType>;

public:
  StencilDataContainerDynImpl() = default;

  StencilDataContainerDynImpl(std::initializer_list<ValueType> il)
    : data_{il}, size_(data_.size())
  {
    for (auto & it : data_){
      ::pressio::ops::set_zero(it);
    }
  }

  StencilDataContainerDynImpl(StencilDataContainerDynImpl const & other) = default;
  StencilDataContainerDynImpl & operator=(StencilDataContainerDynImpl const & other) = default;
  StencilDataContainerDynImpl(StencilDataContainerDynImpl && other) = default;
  StencilDataContainerDynImpl & operator=(StencilDataContainerDynImpl && other) = default;
  ~StencilDataContainerDynImpl() = default;

public:
  std::size_t size() const{ return size_; }

  // non-const overloads
  ValueType & operator()(::pressio::ode::n){
    assert( size_>=1);
    return data_[0];
  }

  ValueType & operator()(::pressio::ode::nMinusOne){
    assert( size_>=2);
    return data_[1];
  }

  ValueType & operator()(::pressio::ode::nMinusTwo){
    assert( size_>=3);
    return data_[2];
  }

  ValueType & operator()(::pressio::ode::nMinusThree){
    assert( size_>=4);
    return data_[3];
  }

  // const overloads
  ValueType const & operator()(::pressio::ode::n) const {
    assert( size_>=1);
    return data_[0];
  }

  ValueType const & operator()(::pressio::ode::nMinusOne) const {
    assert( size_>=2);
    return data_[1];
  }

  ValueType const & operator()(::pressio::ode::nMinusTwo) const {
    assert( size_>=3);
    return data_[2];
  }

  ValueType const & operator()(::pressio::ode::nMinusThree) const {
    assert( size_>=4);
    return data_[3];
  }

private:
  data_type data_;
  std::size_t size_;

};

}}}
#endif  // PRESSIOROM_ODE_IMPL_ODE_STENCIL_DATA_CONTAINER_DYNAMIC_HPP_
