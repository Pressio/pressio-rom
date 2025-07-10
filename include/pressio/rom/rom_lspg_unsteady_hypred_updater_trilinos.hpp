/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_hypred_updater_trilinos.hpp
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

#ifndef PRESSIOROM_ROM_ROM_LSPG_UNSTEADY_HYPRED_UPDATER_TRILINOS_HPP_
#define PRESSIOROM_ROM_ROM_LSPG_UNSTEADY_HYPRED_UPDATER_TRILINOS_HPP_

namespace pressio{ namespace rom{ namespace lspg{

// a = alpha*a + beta*b (a,b potentially non with same distribution)

struct HypRedUpdaterTrilinos
{
  // -----------------
  // TPETRA
  // -----------------
  template<class ScalarType, class ...Args>
  void updateSampleMeshOperandWithStencilMeshOne
  (Tpetra::Vector<Args...> & sample_operand,
   const ScalarType alpha,
   const Tpetra::Vector<Args...> & stencil_operand,
   const ScalarType beta) const
  {

    const auto sample_map   = sample_operand.getMap();
    const auto sample_gIDs  = sample_map->getMyGlobalIndices();
    const auto stencil_map  = stencil_operand.getMap();
    const auto stencil_gIDs = stencil_map->getMyGlobalIndices();

    //loop over LOCAL elements of the sample_operand
    auto sample_op_data  = sample_operand.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    auto stencil_op_data = stencil_operand.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());
    for (std::size_t i=0; i<sample_operand.getLocalLength(); i++)
    {

      // assert that the sample mesh global index is also owned by
      // the stencil map on the calling process
      assert(stencil_map->isNodeGlobalElement(sample_gIDs[i]));

      // we only need to combine things if the global id of the
      // sample operand is found in the stencil so we ask
      // the stencil map what is the local index that corresponds
      // to the global id we are handling
      const auto lid = stencil_map->getLocalElement(sample_gIDs[i]);

      sample_op_data(i,0) = alpha*sample_op_data(i,0)
	+ beta*stencil_op_data(lid,0);
    }
  }

  template<class ScalarType, class ...Args>
  void updateSampleMeshOperandWithStencilMeshOne
  (Tpetra::MultiVector<Args...> & sample_operand,
   const ScalarType alpha,
   const Tpetra::MultiVector<Args...> & stencil_operand,
   const ScalarType beta) const
  {

    const auto sample_map   = sample_operand.getMap();
    const auto sample_gIDs  = sample_map->getMyGlobalIndices();
    const auto stencil_map  = stencil_operand.getMap();
    const auto stencil_gIDs = stencil_map->getMyGlobalIndices();

    //loop over LOCAL elements of the sample_operand
    auto sample_op_data  = sample_operand.getLocalViewHost(Tpetra::Access::ReadWriteStruct());
    auto stencil_op_data = stencil_operand.getLocalViewHost(Tpetra::Access::ReadOnlyStruct());

    for (std::size_t j=0; j<sample_operand.getNumVectors(); j++){
      for (std::size_t i=0; i<sample_operand.getLocalLength(); i++){
	// assert that the sample mesh global index is also owned by
	// the stencil map on the calling process
	assert(stencil_map->isNodeGlobalElement(sample_gIDs[i]));

	// we only need to combine things if the global id of the
	// sample operand is found in the stencil so we ask
	// the stencil map what is the local index that corresponds
	// to the global id we are handling
	const auto lid = stencil_map->getLocalElement(sample_gIDs[i]);

	sample_op_data(i,j) = alpha*sample_op_data(i,j)
	  + beta*stencil_op_data(lid,j);
      }
    }
  }

  // -----------------
  // TPETRA BLOCK
  // -----------------
  template<class ScalarType, class ...Args>
  void updateSampleMeshOperandWithStencilMeshOne
  (Tpetra::BlockVector<Args...> & sample_operand,
   const ScalarType alpha,
   const Tpetra::BlockVector<Args...> & stencil_operand,
   const ScalarType beta) const
  {
    using block_type = Tpetra::BlockVector<Args...>;
    auto sample_op_vv = sample_operand.getVectorView();
    auto stencil_op_vv = const_cast<block_type &>(stencil_operand).getVectorView();
    updateSampleMeshOperandWithStencilMeshOne(sample_op_vv, alpha, stencil_op_vv, beta);
  }

  template<class ScalarType, class ...Args>
  void updateSampleMeshOperandWithStencilMeshOne
  (Tpetra::BlockMultiVector<Args...> & sample_operand,
   const ScalarType alpha,
   const Tpetra::BlockMultiVector<Args...> & stencil_operand,
   const ScalarType beta) const
  {
    auto sample_op_vv = sample_operand.getMultiVectorView();
    auto stencil_op_vv = stencil_operand.getMultiVectorView();
    updateSampleMeshOperandWithStencilMeshOne(sample_op_vv, alpha, stencil_op_vv, beta);
  }

};

}}}
#endif  // PRESSIOROM_ROM_ROM_LSPG_UNSTEADY_HYPRED_UPDATER_TRILINOS_HPP_
