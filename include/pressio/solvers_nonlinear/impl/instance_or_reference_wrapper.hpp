/*
//@HEADER
// ************************************************************************
//
// instance_or_reference_wrapper.hpp
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

#ifndef PRESSIOROM_SOLVERS_NONLINEAR_IMPL_INSTANCE_OR_REFERENCE_WRAPPER_HPP_
#define PRESSIOROM_SOLVERS_NONLINEAR_IMPL_INSTANCE_OR_REFERENCE_WRAPPER_HPP_

namespace pressio{ namespace nlsol{ namespace impl{

template <class T>
class InstanceOrReferenceWrapper;

template <>
class InstanceOrReferenceWrapper<void>
{};

template <class T>
class InstanceOrReferenceWrapper
{
  T value_;

public:
  InstanceOrReferenceWrapper(const InstanceOrReferenceWrapper &) = default;
  InstanceOrReferenceWrapper & operator=(const InstanceOrReferenceWrapper &) = default;
  InstanceOrReferenceWrapper(InstanceOrReferenceWrapper &&) = default;
  InstanceOrReferenceWrapper & operator=(InstanceOrReferenceWrapper &&) = default;
  ~InstanceOrReferenceWrapper() = default;

  template< class _T = T, std::enable_if_t<std::is_default_constructible<_T>::value, int> = 0>
  InstanceOrReferenceWrapper(){}

  explicit InstanceOrReferenceWrapper(T & valIn)
    : value_(valIn){}

  explicit InstanceOrReferenceWrapper(const T & valIn)
    : value_(valIn){}

  explicit InstanceOrReferenceWrapper(T && valIn)
    : value_(std::move(valIn)){}

  T& get() { return value_; }
  T const& get() const { return value_; }
};

template <class T>
class InstanceOrReferenceWrapper<T&>
{
  std::reference_wrapper<T> refObj_;

public:
  InstanceOrReferenceWrapper() = delete;
  InstanceOrReferenceWrapper(const InstanceOrReferenceWrapper &) = default;
  InstanceOrReferenceWrapper & operator=(const InstanceOrReferenceWrapper &) = default;
  InstanceOrReferenceWrapper(InstanceOrReferenceWrapper &&) = default;
  InstanceOrReferenceWrapper & operator=(InstanceOrReferenceWrapper &&) = default;
  ~InstanceOrReferenceWrapper() = default;

  explicit InstanceOrReferenceWrapper(T & valIn) : refObj_(valIn){}

  T& get() { return refObj_.get(); }
  T const& get() const { return refObj_.get(); }
};

template <class T>
class InstanceOrReferenceWrapper<T const &>
{
  std::reference_wrapper<const T> refObj_;

public:
  InstanceOrReferenceWrapper() = delete;
  InstanceOrReferenceWrapper(const InstanceOrReferenceWrapper &) = default;
  InstanceOrReferenceWrapper & operator=(const InstanceOrReferenceWrapper &) = default;
  InstanceOrReferenceWrapper(InstanceOrReferenceWrapper &&) = default;
  InstanceOrReferenceWrapper & operator=(InstanceOrReferenceWrapper &&) = default;
  ~InstanceOrReferenceWrapper() = default;

  explicit InstanceOrReferenceWrapper(const T & valIn) : refObj_(valIn){}

  T const& get() { return refObj_.get(); }
  T const& get() const { return refObj_.get(); }
};

}}} // end of namespace pressio::nlsol::impl
#endif  // PRESSIOROM_SOLVERS_NONLINEAR_IMPL_INSTANCE_OR_REFERENCE_WRAPPER_HPP_
