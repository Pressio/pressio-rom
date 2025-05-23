
#ifndef PRESSIO_ODE_CONCEPTS_ODE_OTHERS_HPP_
#define PRESSIO_ODE_CONCEPTS_ODE_OTHERS_HPP_

namespace pressio{ namespace ode{

template <class T, class = void>
struct Steppable : std::false_type{};

template <class T>
struct Steppable<
  T,
  std::enable_if_t<
    ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval< typename T::state_type & >(),
	std::declval< ::pressio::ode::StepStartAt<typename T::independent_variable_type> >(),
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepSize<typename T::independent_variable_type> >()
	)
       )
      >::value
    //&& impl::stepper_accepting_lvalue_state<T>::value
    >
  > : std::true_type{};

template <class T, class AuxT, class ...Args>
struct SteppableWithAuxiliaryArgs : std::false_type{};

template <class T, class AuxT, class ...Args>
struct SteppableWithAuxiliaryArgs<
  std::enable_if_t<
    ::pressio::has_independent_variable_typedef<T>::value
    && ::pressio::has_state_typedef<T>::value
    && std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval<typename T::state_type & >(),
	std::declval<::pressio::ode::StepStartAt<typename T::independent_variable_type> >(),
	std::declval<::pressio::ode::StepCount >(),
	std::declval<::pressio::ode::StepSize<typename T::independent_variable_type> >(),
	std::declval< AuxT >(), std::declval<Args>()...
	)
       )
      >::value
    //&& impl::variadic_stepper_accepting_lvalue_state<void, T, AuxT, Args...>::value
    >,
  T, AuxT, Args...
  > : std::true_type{};

template <class T>
using StronglySteppable = Steppable<T>;

template <class T, class AuxT, class ...Args>
using StronglySteppableWithAuxiliaryArgs = SteppableWithAuxiliaryArgs<T, AuxT, Args...>;


template <class T, class IndVarType, class StateType, class enable = void>
struct StateObserver : std::false_type{};

template <class T, class IndVarType, class StateType>
struct StateObserver<
  T, IndVarType, StateType,
  std::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T>().operator()
	       (
		std::declval< ::pressio::ode::StepCount >(),
		std::declval< IndVarType >(),
		std::declval<StateType const &>()
		)
	       )
      >::value
    >
  > : std::true_type{};

template <class T, class IndVarType, class Enable = void>
struct StepSizePolicy : std::false_type{};

template <class T, class IndVarType>
struct StepSizePolicy<
  T, IndVarType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepStartAt<IndVarType> >(),
	std::declval< ::pressio::ode::StepSize<IndVarType> & >()
	)
       )
      >::value
    //&& impl::step_size_policy_taking_dt_by_ref<T, IndVarType>::value
    >
  > : std::true_type{};


template <class T, class IndVarType, class Enable = void>
struct StepSizePolicyWithReductionScheme : std::false_type{};

template <class T, class IndVarType>
struct StepSizePolicyWithReductionScheme<
  T, IndVarType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T>()
       (
	std::declval< ::pressio::ode::StepCount >(),
	std::declval< ::pressio::ode::StepStartAt<IndVarType> >(),
	std::declval< ::pressio::ode::StepSize<IndVarType> & >(),
	std::declval< ::pressio::ode::StepSizeMinAllowedValue<IndVarType> & >(),
	std::declval< ::pressio::ode::StepSizeScalingFactor<IndVarType> & >()
	)
       )
      >::value
    //&& impl::step_size_policy_with_reduc_taking_dt_by_ref<T, IndVarType>::value
    >
  > : std::true_type{};

}} // end namespace pressio::ode
#endif  // PRESSIO_ODE_CONCEPTS_ODE_OTHERS_HPP_
