
#ifndef PRESSIO_ROM_LSPG_STEADY_HPP_
#define PRESSIO_ROM_LSPG_STEADY_HPP_

#include "./impl/lspg_steady_system_default.hpp"
#include "./impl/lspg_steady_system_masked.hpp"
#include "./impl/noop.hpp"

namespace pressio{ namespace rom{ namespace lspg{

// -------------------------------------------------------------
// default or hyp-red
// impl-wise, default and hyp-red LSPG are the same
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType>
auto create_steady_problem(const TrialSubspaceType & trialSpace,   /*(1)*/
			   const FomSystemType & fomSystem)

{
  static_assert
    (SteadyFomWithJacobianAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
     "You are trying to create a steady LSPG problem but the \
FOM system does not meet the required SteadyFomWithJacobianAction concept.");

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = impl::NoOperation<void>;
  using system_type = impl::LspgSteadyDefaultSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, scaler_type>;
  return system_type(trialSpace, fomSystem, scaler_type{});
}

// -------------------------------------------------------------
// masked
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType>
auto create_steady_problem(const TrialSubspaceType & trialSpace,  /*(2)*/
			   const FomSystemType & fomSystem,
			   const MaskerType & masker)
{

  static_assert
    (SteadyFomWithJacobianAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
     "You are trying to create a steady LSPG problem but the \
FOM system does not meet the required SteadyFomWithJacobianAction concept.");

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = impl::NoOperation<void>;
  using system_type = impl::LspgSteadyMaskedSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, MaskerType, scaler_type>;
  return system_type(trialSpace, fomSystem, masker, scaler_type{});
}


namespace experimental{
// -------------------------------------------------------------
// default or hyp-red with scaling
// impl-wise, default and hyp-red LSPG are the same
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class ScalingOperatorType>
auto create_steady_problem(const TrialSubspaceType & trialSpace,  /*(3)*/
			   const FomSystemType & fomSystem,
			   const ScalingOperatorType & scaler)
{
  static_assert
    (SteadyFomWithJacobianAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
     "You are trying to create a steady LSPG problem but the \
FOM system does not meet the required SteadyFomWithJacobianAction concept.");

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = std::reference_wrapper<const ScalingOperatorType>;
  using system_type = impl::LspgSteadyDefaultSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, scaler_type>;
  return system_type(trialSpace, fomSystem, scaler);
}

// -------------------------------------------------------------
// masked with scaling
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class ScalingOperatorType>
auto create_steady_problem(const TrialSubspaceType & trialSpace,  /*(4)*/
			   const FomSystemType & fomSystem,
			   const MaskerType & masker,
			   const ScalingOperatorType & scaler)
{
  static_assert
    (SteadyFomWithJacobianAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
     "You are trying to create a steady LSPG problem but the \
FOM system does not meet the required SteadyFomWithJacobianAction concept.");

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = std::reference_wrapper<const ScalingOperatorType>;
  using system_type = impl::LspgSteadyMaskedSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, MaskerType, scaler_type>;
  return system_type(trialSpace, fomSystem, masker, scaler);
}

} // end experimental

}}} // end pressio::rom::lspg
#endif  // PRESSIO_ROM_LSPG_STEADY_HPP_
