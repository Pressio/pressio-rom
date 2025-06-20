
#ifndef PRESSIO_ROM_IMPL_GALERKIN_HELPERS_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_HELPERS_HPP_

namespace pressio{ namespace rom{ namespace impl{

template<class T = void>
void valid_scheme_for_explicit_galerkin_else_throw(::pressio::ode::StepScheme name,
						   const std::string & str){
  if (!::pressio::ode::is_explicit_scheme(name)){
    throw std::runtime_error(str + " requires an explicit stepper");
  }
}

template<class T = void>
void valid_scheme_for_implicit_galerkin_else_throw(::pressio::ode::StepScheme name,
						   const std::string & str){
  if (!::pressio::ode::is_implicit_scheme(name)){
    throw std::runtime_error(str + " requires an implicit stepper");
  }
}

// --------------------------------------------------------------
// CreateGalerkinRhs
// --------------------------------------------------------------
template<class T, class = void> struct CreateGalerkinRhs;

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct CreateGalerkinRhs<
  T, std::enable_if_t< ::pressio::is_vector_eigen<T>::value >
  >
{
  T operator()(std::size_t ext){ return T(ext); }
};
#endif

// --------------------------------------------------------------
// CreateGalerkinMassMatrix
// --------------------------------------------------------------
template<class T, class = void>
struct CreateGalerkinMassMatrix;

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct CreateGalerkinMassMatrix<
  T, std::enable_if_t< ::pressio::is_dense_matrix_eigen<T>::value >
  >
{
  T operator()(std::size_t ext){ return T(ext, ext); }
};
#endif

// ------------------------------------------
// this is an alias becuase it can done in the same way
template<class JacType, class = void>
using CreateGalerkinJacobian = CreateGalerkinMassMatrix<JacType>;


template<class TrialSubspaceType, class FomSystemType>
void galerkin_static_check_trial_and_system(const TrialSubspaceType & trialSpace,
				     const FomSystemType & fomSystem)
{
  static_assert(PossiblyAffineTrialColumnSubspace<TrialSubspaceType>::value,
     "You are trying to create a steady galerkin problem but the \
trialSpace does not meet the required PossiblyAffineTrialColumnSubspace concept.");

  static_assert
    (SteadyFomWithJacobianAction<
     FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
     "You are trying to create a steady galerkin problem but the \
FOM system does not meet the required SteadyFomWithJacobianAction concept.");
}

template<typename T>
void steady_galerkin_static_check_api_return_type(){
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
#endif  // PRESSIO_ROM_IMPL_GALERKIN_HELPERS_HPP_
