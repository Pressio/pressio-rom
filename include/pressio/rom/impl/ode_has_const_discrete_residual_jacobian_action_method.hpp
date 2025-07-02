
#ifndef PRESSIOROM_ROM_IMPL_ODE_HAS_CONST_DISCRETE_RESIDUAL_JACOBIAN_ACTION_METHOD_HPP_
#define PRESSIOROM_ROM_IMPL_ODE_HAS_CONST_DISCRETE_RESIDUAL_JACOBIAN_ACTION_METHOD_HPP_

namespace pressio{ namespace rom{

template <
  class T,
  int n,
  class StepType,
  class IndVarType,
  class StateType,
  class ResidualType,
  class ManifoldJacobian,
  class JacobianType,
  class = void
  >
struct has_const_discrete_residual_jacobian_action_method : std::false_type{};

template <
  class T, class StepType, class IndVarType, class StateType,
  class ResidualType, class ManifoldJacobian, class JacobianType>
struct has_const_discrete_residual_jacobian_action_method<
  T, 1, StepType, IndVarType, StateType, ResidualType, ManifoldJacobian, JacobianType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const>().discreteTimeResidualAndJacobianAction
       (
	  std::declval<StepType const &>(),
	  std::declval<IndVarType const &>(),
	  std::declval<IndVarType const &>(),
	  std::declval<ResidualType &>(),
	  std::declval<ManifoldJacobian const &>(),
    std::declval< std::optional<JacobianType*> >(),
	  std::declval<StateType const&>()
	)
       )
      >::value
    >
  > : std::true_type{};

template <
  class T, class StepType, class IndVarType, class StateType,
  class ResidualType, class ManifoldJacobian, class JacobianType>
struct has_const_discrete_residual_jacobian_action_method<
  T, 2, StepType, IndVarType, StateType, ResidualType, ManifoldJacobian, JacobianType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
      std::declval<T const>().discreteTimeResidualAndJacobianAction
      (
       std::declval<StepType const &>(),
       std::declval<IndVarType const &>(),
       std::declval<IndVarType const &>(),
       std::declval<ResidualType &>(),
       std::declval<ManifoldJacobian const &>(),
       std::declval< std::optional<JacobianType*> >(),
       std::declval<StateType const&>(),
       std::declval<StateType const&>()
       )
       )
      >::value
  >
  > : std::true_type{};

template <
  class T, class StepType, class IndVarType, class StateType,
  class ResidualType, class ManifoldJacobian, class JacobianType>
struct has_const_discrete_residual_jacobian_action_method<
  T, 3, StepType, IndVarType, StateType, ResidualType, ManifoldJacobian, JacobianType,
  std::enable_if_t<
    std::is_void<
      decltype
      (
      std::declval<T const>().discreteTimeResidualAndJacobianAction
      (
       std::declval<StepType const &>(),
       std::declval<IndVarType const &>(),
       std::declval<IndVarType const &>(),
       std::declval<ResidualType &>(),
       std::declval<ManifoldJacobian const &>(),
       std::declval< std::optional<JacobianType*> >(),
       std::declval<StateType const&>(),
       std::declval<StateType const&>(),
       std::declval<StateType const&>()
       )
       )
  >::value
  >
  > : std::true_type{};

}}
#endif  // PRESSIOROM_ROM_IMPL_ODE_HAS_CONST_DISCRETE_RESIDUAL_JACOBIAN_ACTION_METHOD_HPP_
