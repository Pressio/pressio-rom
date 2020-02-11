/*
//@HEADER
// ************************************************************************
//
// containers_vector_sharedmem_blaze_dynamic.hpp
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

#ifdef PRESSIO_ENABLE_TPL_BLAZE
#ifndef CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_BLAZE_DYNAMIC_HPP_
#define CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_BLAZE_DYNAMIC_HPP_

#include "../../base/containers_container_base.hpp"
#include "../../base/containers_container_sharedmem_base.hpp"
#include "../../base/containers_vector_sharedmem_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<wrapped_type,
    ::pressio::mpl::enable_if_t<
      containers::meta::is_dynamic_vector_blaze<wrapped_type>::value
      >
    >
  : public ContainerBase< Vector<wrapped_type>, wrapped_type >,
    public ContainerSharedMemBase< Vector<wrapped_type> >,
    public VectorSharedMemBase< Vector<wrapped_type> >
{

  using this_t = Vector<wrapped_type>;
  using mytraits = typename details::traits<this_t>;
  using sc_t = typename mytraits::scalar_t;
  using ord_t = typename  mytraits::ordinal_t;
  using wrap_t = typename mytraits::wrapped_t;

public:
  Vector() = default;

  explicit Vector(ord_t insize){
    this->resize(insize);
  }
  explicit Vector(const wrap_t & src) : data_(src){}

public:
  sc_t & operator [] (ord_t i){
    return data_[i];
  };

  sc_t const & operator [] (ord_t i) const{
    return data_[i];
  };

  sc_t & operator()(ord_t i){
    return data_[i];
  };
  sc_t const & operator()(ord_t i) const{
    return data_[i];
  };

  // compound assignment when type(b) = type(this)
  // this += b
  // template <typename T,
  // 	    ::pressio::mpl::enable_if_t<
  // 	      std::is_same<T,this_t>::value> * = nullptr>
  this_t & operator+=(const this_t & other) {
    assert( other.size() == this->size() );
    this->data_ += *other.data();
    return *this;
  }

  // compound assignment when type(b) = type(this)
  // this -= b
  // template <typename T,
  // 	    ::pressio::mpl::enable_if_t<
  // 	      std::is_same<T,this_t>::value> * = nullptr>
  this_t & operator-=(const this_t & other) {
    assert( other.size() == this->size() );
    this->data_ -= *other.data();
    return *this;
  }

private:
  wrap_t const * dataImpl() const{
    return &data_;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  bool emptyImpl() const{
    return this->size()==0 ? true : false;
  }

  ord_t extentImpl() const {
    return data_.size();
  }

private:
  friend ContainerBase< this_t, wrapped_type >;
  friend ContainerSharedMemBase< this_t >;
  friend VectorSharedMemBase< this_t >;

private:
  wrap_t data_ = {};

};//end class
}}//end namespace pressio::containers

#endif /* CONTAINERS_VECTOR_CONCRETE_VECTOR_SHAREDMEM_BLAZE_DYNAMIC_HPP_ */
#endif //PRESSIO_ENABLE_TPL_BLAZE
