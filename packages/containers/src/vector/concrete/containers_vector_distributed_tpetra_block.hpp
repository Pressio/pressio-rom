/*
//@HEADER
// ************************************************************************
//
// containers_vector_distributed_tpetra_block.hpp
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

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#ifndef CONTAINERS_VECTOR_CONCRETE_VECTOR_DISTRIBUTED_TPETRA_BLOCK_HPP_
#define CONTAINERS_VECTOR_CONCRETE_VECTOR_DISTRIBUTED_TPETRA_BLOCK_HPP_

#include "../../base/containers_container_base.hpp"
#include "../../base/containers_container_distributed_base.hpp"
#include "../../base/containers_vector_distributed_base.hpp"

namespace pressio{ namespace containers{

template <typename wrapped_type>
class Vector<
  wrapped_type,
  ::pressio::mpl::enable_if_t<
    meta::is_vector_tpetra_block<
      wrapped_type
      >::value
    >
  >
  : public ContainerBase< Vector<wrapped_type>, wrapped_type >,
    public ContainerDistributedBase< Vector<wrapped_type> >,
    public VectorDistributedBase< Vector<wrapped_type> >
{

  using this_t		= Vector<wrapped_type>;
  using der_t		= this_t;
  using sc_t		= typename details::traits<this_t>::scalar_t;
  using LO_t		= typename details::traits<this_t>::local_ordinal_t;
  using GO_t		= typename details::traits<this_t>::global_ordinal_t;
  using device_t	= typename details::traits<this_t>::device_t;
  using wrap_t		= typename details::traits<this_t>::wrapped_t;
  using map_t		= typename details::traits<this_t>::data_map_t;
  using mpicomm_t	= typename details::traits<this_t>::communicator_t;

public:
  Vector() = delete;

  /* Block MV/V still missing a copy constructor,
   * see https://github.com/trilinos/Trilinos/issues/4627
   * so for now we construct this one using other's data */

  explicit Vector(const wrap_t & vecobj)
    : data_( *vecobj.getMap(),
  	     vecobj.getBlockSize()){
    // just a trick to copy data
    data_.update(::pressio::utils::constants::one<sc_t>(),
		 vecobj,
		 ::pressio::utils::constants::zero<sc_t>());
  }

  // here we do not default the copy and move because if we did that,
  // it would use the tpetra copy/move which have view semantics
  // which is not what we want here (for the time being)

  // copy cnstr delegating (for now) to the one above
  Vector(Vector const & other) : Vector(*other.data()){}

  // copy assignment
  Vector & operator=(const Vector & other){
    if (&other != this){
      this->data_.update(::pressio::utils::constants::one<sc_t>(),
			 *other.data(),
			 ::pressio::utils::constants::zero<sc_t>() );
    }
    return *this;
  }

  // move cnstr
  Vector(Vector && other) : Vector(*other.data()){}

  // move assignment
  Vector & operator=(Vector && other){
    this->data_.update(::pressio::utils::constants::one<sc_t>(),
		       *other.data(),
		       ::pressio::utils::constants::zero<sc_t>() );
    return *this;
  }

  ~Vector() = default;

public:
  // compound add assignment when type(b) = type(this)
  // this += b
  this_t & operator+=(const this_t & other) {
    this->data_.update(::pressio::utils::constants::one<sc_t>(),
		       *other.data(),
		       ::pressio::utils::constants::one<sc_t>() );
    return *this;
  }

  // compound add assignment when type(b) = type(this)
  // this -= b
  this_t & operator-=(const this_t & other) {
    this->data_.update(::pressio::utils::constants::negOne<sc_t>(),
		       *other.data(),
		       ::pressio::utils::constants::one<sc_t>() );
    return *this;
  }

private:

  // map_t const & getDataMapImpl() const{
  //   return *data_.getMap();
  // }

  wrap_t const * dataImpl() const{
    return &data_;
  }

  wrap_t * dataImpl(){
    return &data_;
  }

  wrap_t dataCpImpl(){
    return data_;
  }

  GO_t extentImpl(std::size_t i) const{
    assert(i==0);
    return data_.getMap()->getGlobalNumElements();
  }

  LO_t extentLocalImpl(std::size_t i) const{
    assert(i==0);
    return data_.getMap()->getNodeNumElements();
  }

private:
  friend ContainerBase< this_t, wrapped_type >;
  friend ContainerDistributedBase< this_t >;
  friend VectorDistributedBase< this_t >;

private:
  wrap_t data_ = {};

};//end class

}}//end namespace pressio::containers
#endif
#endif
