
#ifndef PRESSIOROM_ROM_IMPL_LSPG_HELPERS_HPP_
#define PRESSIOROM_ROM_IMPL_LSPG_HELPERS_HPP_

namespace pressio{ namespace rom{ namespace impl{

template<class T= void>
void valid_scheme_for_lspg_else_throw(::pressio::ode::StepScheme name){
  if (   name != ::pressio::ode::StepScheme::BDF1
      && name != ::pressio::ode::StepScheme::BDF2)
  {
    throw std::runtime_error("LSPG currently accepting BDF1 or BDF2");
  }
}

template<class TrialSubspaceType, class FomSystemType>
void lspg_static_check_trial_and_system(const TrialSubspaceType & trialSpace,
					const FomSystemType & fomSystem)
{
  static_assert(PossiblyAffineTrialColumnSubspace<TrialSubspaceType>::value,
     "You are trying to create a steady lspg problem but the \
trialSpace does not meet the required PossiblyAffineTrialColumnSubspace concept.");

  static_assert
    (SteadyFomWithJacobianAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
     "You are trying to create a steady lspg problem but the \
FOM system does not meet the required SteadyFomWithJacobianAction concept.");
}

template<typename T>
void steady_lspg_static_check_api_return_type(){
  static constexpr bool val =
#ifdef PRESSIO_ENABLE_CXX20
    nonlinearsolvers::NonlinearSystemFusingResidualAndJacobian<T>;
#else
  nonlinearsolvers::NonlinearSystemFusingResidualAndJacobian<T>::value;
#endif
  static_assert(val,
		"The return type must satisify the NonlinearSystemFusingResidualAndJacobian concept.");
}

}}} // end pressio::rom::impl
#endif  // PRESSIOROM_ROM_IMPL_LSPG_HELPERS_HPP_
