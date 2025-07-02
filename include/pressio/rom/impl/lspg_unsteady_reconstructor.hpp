/*
//@HEADER
// ************************************************************************
//
// lspg_unsteady_reconstructor.hpp
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

#ifndef PRESSIOROM_ROM_IMPL_LSPG_UNSTEADY_RECONSTRUCTOR_HPP_
#define PRESSIOROM_ROM_IMPL_LSPG_UNSTEADY_RECONSTRUCTOR_HPP_

#include <fstream>

namespace pressio{ namespace rom{ namespace impl{

template<typename reduced_state_type>
auto read_rom_states_and_times_from_ascii(std::string const & filein, std::size_t ext)
{
  using scalar_type = scalar_trait_t<reduced_state_type>;

  std::ifstream reduced_state_file(filein);
  assert(reduced_state_file.good());

  std::string line;
  std::vector<scalar_type> times;
  std::vector<reduced_state_type> reduced_states;
  while (std::getline(reduced_state_file, line)){
    std::istringstream stm(line);
    scalar_type time;
    reduced_state_type reduced_state(ext);
    scalar_type reduced_state_i;
    stm >> time;
    times.push_back(time);

    int c=0;
    while (stm >> reduced_state_i){
      reduced_state[c++] = reduced_state_i;
    }
    reduced_states.push_back(reduced_state);
  }
  reduced_state_file.close();
  return std::make_tuple(times, reduced_states);
}

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template<typename ViewT>
auto _rank1_view_to_stdvector(ViewT view)
{
  using sc_t = std::remove_cv_t<scalar_trait_t<ViewT>>;
  std::vector<sc_t> v(view.extent(0));
  for (std::size_t k=0; k<view.extent(0); ++k){
    v[k] = view(k,0);
  }
  return v;
}

template<typename ViewT>
auto _rank2_view_to_stdvector(ViewT view)
{
  using sc_t = std::remove_cv_t<scalar_trait_t<ViewT>>;
  std::vector<sc_t> v(view.extent(0)*view.extent(1));
  std::size_t k=0;
  for (std::size_t j=0; j<view.extent(1); ++j){
    for (std::size_t i=0; i<view.extent(0); ++i){
      v[k++] = view(i,j);
    }
  }
  return v;
}
#endif

template<typename sc_t, typename size_t>
void write_vector_to_binary(const std::string filename,
			    const sc_t * A,
			    size_t n)
{
  std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
  out.write((char*) A, n*sizeof(sc_t) );
  out.close();
}

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
struct DefaultWriter
{
  std::string prependToFileOut_;

  template<class Rt, class Jt>
  void operator()(std::size_t i, const Rt & R, const Jt & JPhi) const
  {
    const auto myRank = R.getMap()->getComm()->getRank();
    const std::string finalPart = "rank_" + std::to_string(myRank) + "_step_" + std::to_string(i) + ".bin";

    auto r_view = R.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto r_stdv = _rank1_view_to_stdvector(r_view);
    const std::string R_f = prependToFileOut_ + "residual_" + finalPart;
    write_vector_to_binary(R_f, r_stdv.data(), r_stdv.size());

    auto jphi_view = JPhi.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto jphi_stdv = _rank2_view_to_stdvector(jphi_view);
    std::string Jphi_f = prependToFileOut_ + "jacobian_action_" + finalPart;
    write_vector_to_binary(Jphi_f, jphi_stdv.data(), jphi_stdv.size());
  }
};

template<class MapType>
void write_map_to_file(MapType const & map)
{
  std::ofstream outMapFile("row_map.txt");
  auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(outMapFile));
  map.describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
  outMapFile.close();
}

template<class T, typename ResidualType, typename TrialSubspaceType, typename enable = void>
struct is_writer : std::false_type {};

template<typename T, typename ResidualType, typename TrialSubspaceType>
struct is_writer<
  T, ResidualType, TrialSubspaceType,
  std::enable_if_t<
    std::is_void<
     decltype(
       std::declval<T const>()(
         std::size_t{},
         std::declval<ResidualType const>(),
         std::declval<typename TrialSubspaceType::basis_matrix_type const>()
        )
     )
    >::value
  >
> : std::true_type{};

#endif

template <class TrialSubspaceType>
class LspgReconstructor{
  using rom_state_type = typename TrialSubspaceType::reduced_state_type;
  using fom_state_type = typename TrialSubspaceType::full_state_type;
  using scalar_type = scalar_trait_t<rom_state_type>;

  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;

public:
  LspgReconstructor(const TrialSubspaceType & trialSubspace)
    : trialSubspace_(trialSubspace){}

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  //
  // constrained for fully discrete system
  //
  template <
    std::size_t n,
    class FomSystemType,
    class _TrialSubspaceType = TrialSubspaceType
    >
  std::enable_if_t<
    RealValuedFullyDiscreteSystemWithJacobianAction<
      FomSystemType, n, typename _TrialSubspaceType::basis_matrix_type>::value
    >
  execute(
    const FomSystemType & fomSystem,
	  const std::string & romStateFilename,
	  std::optional<std::string> filenamePrepend = {}) const
  {
    static_assert(n==2,
    "lspg reconstructor for a fully discrete system currently supports TotalNumberOfDesiredStates==2");

    DefaultWriter writer{ filenamePrepend.value_or("") };
    execute_for_fully_discrete_time_impl(fomSystem, romStateFilename, writer);
  }

  template <
    std::size_t n,
    class FomSystemType,
    class CustomWriterType,
    class _TrialSubspaceType = TrialSubspaceType
    >
  std::enable_if_t<
    RealValuedFullyDiscreteSystemWithJacobianAction<
      FomSystemType, n, typename _TrialSubspaceType::basis_matrix_type>::value
      && is_writer<CustomWriterType, typename FomSystemType::discrete_residual_type, _TrialSubspaceType>::value
    >
  execute(
    const FomSystemType & fomSystem,
	  const std::string & romStateFilename,
    const CustomWriterType & writer) const
  {
    static_assert(n==2,
    "lspg reconstructor for a fully discrete system currently supports TotalNumberOfDesiredStates==2");

    execute_for_fully_discrete_time_impl(fomSystem, romStateFilename, writer);
  }

  //
  // constrained for semi-discrete system
  //
  template <
    class FomSystemType,
    class _TrialSubspaceType = TrialSubspaceType
    >
  std::enable_if_t<
    RealValuedSemiDiscreteFomWithJacobianAction<
      FomSystemType, typename _TrialSubspaceType::basis_matrix_type
      >::value
    >
  execute(
    const FomSystemType & fomSystem,
	  std::string const & romStateFilename,
	  ::pressio::ode::StepScheme schemeName,
	  std::optional<std::string> filenamePrepend = {}) const
  {
    DefaultWriter writer{ filenamePrepend.value_or("") };
    execute_for_semi_discrete_time_impl(fomSystem, romStateFilename, schemeName, writer);
  }

  template <
    class FomSystemType,
    class CustomWriterType,
    class _TrialSubspaceType = TrialSubspaceType
    >
  std::enable_if_t<
    RealValuedSemiDiscreteFomWithJacobianAction<
      FomSystemType, typename _TrialSubspaceType::basis_matrix_type>::value
      && is_writer<CustomWriterType, typename FomSystemType::rhs_type, _TrialSubspaceType>::value
    >
  execute(
    const FomSystemType & fomSystem,
	  std::string const & romStateFilename,
	  ::pressio::ode::StepScheme schemeName,
    const CustomWriterType & writer) const
  {
    execute_for_semi_discrete_time_impl(fomSystem, romStateFilename, schemeName, writer);
  }


private:
  template <class FomSystemType, class WriterType>
  void execute_for_fully_discrete_time_impl(
      const FomSystemType & fomSystem,
      const std::string & romStateFilename,
      const WriterType & writer) const
  {
    const auto & trialSub = trialSubspace_.get();
    const auto & phi = trialSub.basisOfTranslatedSpace();

    // 1. allocate what we need
    auto state_np1 = trialSub.createFullState();
    auto state_n   = trialSub.createFullState();
    auto R = fomSystem.createDiscreteTimeResidual();
    auto JTimesPhi = fomSystem.createResultOfDiscreteTimeJacobianActionOn(phi);
    assert( *R.getMap() == *JTimesPhi.getMap() );

    // 2. write the row map
    write_map_to_file(*R.getMap());

    // 3. read states
    const std::size_t numModes = trialSubspace_.get().dimension();
    auto [times, reducedStates] = read_rom_states_and_times_from_ascii<rom_state_type>(romStateFilename, numModes);

    trialSub.mapFromReducedState(reducedStates[0], state_n);
    for (std::size_t i = 1; i < times.size(); i++){
      const auto t_n   = times[i-1];
      const auto t_np1 = times[i];
      const auto dt    = t_np1 - t_n;

      trialSub.mapFromReducedState(reducedStates[i], state_np1);
      fomSystem.discreteTimeResidualAndJacobianAction(i, t_np1, dt, R,
						      phi, &JTimesPhi, state_np1, state_n);
      writer(i, R, JTimesPhi);
      pressio::ops::deep_copy(state_n, state_np1);
    }
  }

  template <class FomSystemType, class WriterType>
  void execute_for_semi_discrete_time_impl(
      const FomSystemType & fomSystem,
      const std::string & romStateFilename,
      ::pressio::ode::StepScheme schemeName,
      const WriterType & writer) const
  {
    assert(schemeName == pressio::ode::StepScheme::BDF1);

    const auto & trialSub = trialSubspace_.get();
    const auto & phi = trialSub.basisOfTranslatedSpace();

    // 1. allocate what we need
    auto state_np1 = trialSub.createFullState();
    // for bdf1, this contains state_n
    ode::ImplicitStencilStatesDynamicContainer<fom_state_type> fomStencilStates{trialSub.createFullState()};

    auto R = fomSystem.createRhs();
    auto JTimesPhi = fomSystem.createResultOfJacobianActionOn(phi);
    assert( *R.getMap() == *JTimesPhi.getMap() );

    // 2. write out the row map
    write_map_to_file(*R.getMap());

    // 3. read states
    const std::size_t numModes = trialSubspace_.get().dimension();
    auto [times, reducedStates] = read_rom_states_and_times_from_ascii<rom_state_type>(romStateFilename, numModes);

    auto & state_n = fomStencilStates(::pressio::ode::n());
    const auto one = static_cast<scalar_type>(1);
    trialSub.mapFromReducedState(reducedStates[0], state_n);
    for (std::size_t i = 1; i < times.size(); i++){
      const auto t_n   = times[i-1];
      const auto t_np1 = times[i];
      const auto dt    = t_np1 - t_n;

      trialSub.mapFromReducedState(reducedStates[i], state_np1);
      fomSystem.rhs(state_np1, t_np1, R);
      fomSystem.applyJacobian(state_np1, phi, t_np1, JTimesPhi);

      ::pressio::ode::impl::discrete_residual(pressio::ode::BDF1(), state_np1, R, fomStencilStates, dt);
      const auto factor = dt*::pressio::ode::constants::bdf1<scalar_type>::c_f_;
      ::pressio::ops::update(JTimesPhi, factor, phi, one);

      writer(i, R, JTimesPhi);
      pressio::ops::deep_copy(state_n, state_np1);
    }
  }
#endif

};

}}}
#endif  // PRESSIOROM_ROM_IMPL_LSPG_UNSTEADY_RECONSTRUCTOR_HPP_
