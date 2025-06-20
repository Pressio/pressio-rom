
#ifndef PRESSIO_ROM_GALERKIN_STEADY_HPP_
#define PRESSIO_ROM_GALERKIN_STEADY_HPP_

#include "./impl/galerkin_helpers.hpp"
#include "./impl/galerkin_steady_system_default.hpp"
#include "./impl/galerkin_steady_system_hypred.hpp"
#include "./impl/galerkin_steady_system_masked.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

// ------------------------------------------------------------------------
// default
// ------------------------------------------------------------------------

template<
  class TrialSubspaceType,
  class FomSystemType>
auto create_steady_problem(const TrialSubspaceType & trialSpace,   /*(1)*/
			   const FomSystemType & fomSystem)
{

  static_assert
    (SteadyFomWithJacobianAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
     "You are trying to create a steady galerkin problem but the \
FOM system does not meet the required SteadyFomWithJacobianAction concept.");

  using reduced_state_t    = typename TrialSubspaceType::reduced_state_type;
  using reduced_r_t = impl::steady_galerkin_default_reduced_residual_t<TrialSubspaceType>;
  using reduced_j_t = impl::steady_galerkin_default_reduced_jacobian_t<TrialSubspaceType>;
  using return_type = impl::GalerkinSteadyDefaultSystem<
    reduced_state_t, reduced_r_t, reduced_j_t,
    TrialSubspaceType, FomSystemType>;
  return return_type(trialSpace, fomSystem);
}

// ------------------------------------------------------------------------
// hyper-reduced
// ------------------------------------------------------------------------

template<
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType>
auto create_steady_problem(const TrialSubspaceType & trialSpace,   /*(2)*/
			   const FomSystemType & fomSystem,
			   const HyperReducerType & hyperReducer)
{

  static_assert
    (SteadyFomWithJacobianAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
     "You are trying to create a steady galerkin problem but the \
FOM system does not meet the required SteadyFomWithJacobianAction concept.");

  using reduced_state_t = typename TrialSubspaceType::reduced_state_type;
  using reduced_r_t = impl::steady_galerkin_default_reduced_residual_t<TrialSubspaceType>;
  using reduced_j_t = impl::steady_galerkin_default_reduced_jacobian_t<TrialSubspaceType>;
  using return_type = impl::GalerkinSteadyHypRedSystem<
    reduced_state_t, reduced_r_t, reduced_j_t,
    TrialSubspaceType, FomSystemType, HyperReducerType>;
  return return_type(trialSpace, fomSystem, hyperReducer);
}

// ------------------------------------------------------------------------
// hyper-reduced masked
// ------------------------------------------------------------------------

template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType>
auto create_steady_problem(const TrialSubspaceType & trialSpace,   /*(3)*/
			   const FomSystemType & fomSystem,
			   const MaskerType & masker,
			   const HyperReducerType & hyperReducer)
{

  static_assert
    (SteadyFomWithJacobianAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
     "You are trying to create a steady galerkin problem but the \
FOM system does not meet the required SteadyFomWithJacobianAction concept.");

  using reduced_state_t = typename TrialSubspaceType::reduced_state_type;
  using reduced_r_t = impl::steady_galerkin_default_reduced_residual_t<TrialSubspaceType>;
  using reduced_j_t = impl::steady_galerkin_default_reduced_jacobian_t<TrialSubspaceType>;
  using return_type = impl::GalerkinSteadyMaskedSystem<
    reduced_state_t, reduced_r_t, reduced_j_t,
    TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>;
  return return_type(trialSpace, fomSystem, masker, hyperReducer);
}

}}} // end pressio::rom::galerkin
#endif  // PRESSIO_ROM_GALERKIN_STEADY_HPP_
