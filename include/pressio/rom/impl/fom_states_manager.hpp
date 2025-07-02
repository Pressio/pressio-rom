/*
//@HEADER
// ************************************************************************
//
// rom_manager_fom_states.hpp
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

#ifndef PRESSIOROM_ROM_IMPL_FOM_STATES_MANAGER_HPP_
#define PRESSIOROM_ROM_IMPL_FOM_STATES_MANAGER_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <class TrialSubspaceType>
class FomStatesManager
{
  using fom_state_type = typename TrialSubspaceType::full_state_type;

public:
  using data_type  = std::vector<fom_state_type>;
  using value_type = fom_state_type;

  FomStatesManager() = delete;
  FomStatesManager(const FomStatesManager &) = default;
  FomStatesManager & operator=(const FomStatesManager &) = delete;
  FomStatesManager(FomStatesManager &&) = default;
  FomStatesManager & operator=(FomStatesManager &&) = delete;
  ~FomStatesManager() = default;

  FomStatesManager(const TrialSubspaceType & trialSubspace,
		       std::initializer_list<fom_state_type> il)
    : trialSubspace_(trialSubspace), data_(il)
  {
    this->setZero();
  }

public:
  std::size_t size() const { return data_.size(); }

  // n+1
  fom_state_type const & operator()(::pressio::ode::nPlusOne) const {
    assert(data_.size() >=1); return data_[0];
  }

  // n
  fom_state_type const & operator()(::pressio::ode::n) const {
    assert(data_.size() >=2); return data_[1];
  }

  // n-1
  fom_state_type const & operator()(::pressio::ode::nMinusOne) const {
    assert(data_.size() >=3); return data_[2];
  }

  // n-2
  fom_state_type const & operator()(::pressio::ode::nMinusTwo) const {
    assert(data_.size() >=4); return data_[3];
  }

  // n+1
  template <class RomStateType>
  void reconstructAtWithoutStencilUpdate(const RomStateType & romStateIn,
				    ::pressio::ode::nPlusOne /*tag*/){
    assert(data_.size() >=1);
    trialSubspace_.get().mapFromReducedState(romStateIn, data_[0]);
  }

  // n
  template <class RomStateType>
  void reconstructAtWithoutStencilUpdate(const RomStateType & romStateIn,
				    ::pressio::ode::n /*tag*/){
    assert(data_.size() >=2);
    trialSubspace_.get().mapFromReducedState(romStateIn, data_[1]);
  }

  // n-1
  template <class RomStateType>
  void reconstructAtWithoutStencilUpdate(const RomStateType & romStateIn,
				    ::pressio::ode::nMinusOne /*tag*/){
    assert(data_.size() >=3);
    trialSubspace_.get().mapFromReducedState(romStateIn, data_[2]);
  }

  // n-2
  template <class RomStateType>
  void reconstructAtWithoutStencilUpdate(const RomStateType & romStateIn,
				    ::pressio::ode::nMinusTwo /*tag*/){
    assert(data_.size() >=4);
    trialSubspace_.get().mapFromReducedState(romStateIn, data_[4]);
  }

  template <class RomStateType>
  void reconstructAtWithStencilUpdate(const RomStateType & romStateIn,
				      ::pressio::ode::n /*tag*/)
  {

    assert(data_.size() >=2);

    if (data_.size() == 2){
      /* when n == 2, it means I only have n+1 and n
       * so reconstructing at n, we just overwrite data_[1] */
      trialSubspace_.get().mapFromReducedState(romStateIn, data_[1]);
    }
    else if (data_.size() == 3){
      /* when n == 3, we have y_n+1, y_n, y_n-1 */
      // y_n becomes  y_n-1
      ::pressio::ops::deep_copy(data_[2], data_[1]);
      // reconstruct y_n
      trialSubspace_.get().mapFromReducedState(romStateIn, data_[1]);
    }
    else if (data_.size() == 4){
      /* when n == 4, we have y_n+1, y_n, y_n-1, y_n-2 */

      // y_n-1 becomes y_n-2
      ::pressio::ops::deep_copy(data_[3], data_[2]);

      // y_n   becomes y_n-1
      ::pressio::ops::deep_copy(data_[2], data_[1]);

      // reconstruct y_n
      trialSubspace_.get().mapFromReducedState(romStateIn, data_[1]);
    }
  }

private:
  void setZero(){
    for (std::size_t i=0; i<data_.size(); i++)
      ::pressio::ops::set_zero(data_[i]);
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  data_type data_;
};

}}}//end namespace pressio::rom::impl

#endif  // PRESSIOROM_ROM_IMPL_FOM_STATES_MANAGER_HPP_
