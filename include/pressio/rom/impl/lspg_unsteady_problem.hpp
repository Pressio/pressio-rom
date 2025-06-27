
#ifndef PRESSIO_ROM_IMPL_LSPG_UNSTEADY_PROBLEM_HPP_
#define PRESSIO_ROM_IMPL_LSPG_UNSTEADY_PROBLEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <int TotalNumberOfDesiredStates, class = void>
struct DeducedStepperType;

template <> struct DeducedStepperType<-1>{
  // note: to deduce the stepper_type it does not really matter
  // what scheme enum value we use, as long as it is an implicit one
  template<class T>
  using type = decltype(::pressio::ode::create_implicit_stepper_with_custom_policy
			(::pressio::ode::StepScheme::BDF1,
			 std::declval<T &>()
			 ));
};

template <int TotalNumberOfDesiredStates>
struct DeducedStepperType<
  TotalNumberOfDesiredStates,
  std::enable_if_t< (TotalNumberOfDesiredStates>0) >
  >
{
  template<class T>
  using type = decltype(::pressio::ode::create_implicit_stepper<
			TotalNumberOfDesiredStates>(std::declval<T &>()));
};

template <
  int TotalNumberOfDesiredStates,
  class TrialSubspaceType,
  class ResJacPolicyOrFullyDiscreteSystemType
  >
class LspgUnsteadyProblem
{

  using stepper_type = typename DeducedStepperType<TotalNumberOfDesiredStates>::template
    type<ResJacPolicyOrFullyDiscreteSystemType>;

  using fom_states_manager_type = FomStatesManager<TrialSubspaceType>;

public:
  using independent_variable_type  = typename ResJacPolicyOrFullyDiscreteSystemType::independent_variable_type;
  using state_type    = typename TrialSubspaceType::reduced_state_type;
  using residual_type = typename ResJacPolicyOrFullyDiscreteSystemType::residual_type;
  using jacobian_type = typename ResJacPolicyOrFullyDiscreteSystemType::jacobian_type;

  template<class FomSystemType, class ...Args>
  LspgUnsteadyProblem(const TrialSubspaceType & trialSubspace,
		      const FomSystemType & fomSystem,
		      Args && ... args)
    : fomStatesManager_(std::make_unique<FomStatesManager<TrialSubspaceType>>(create_lspg_fom_states_manager<TotalNumberOfDesiredStates>(trialSubspace))),
      rjPolicyOrFullyDiscreteSystem_(std::make_unique<ResJacPolicyOrFullyDiscreteSystemType>(trialSubspace, fomSystem,
											     *fomStatesManager_,
											     std::forward<Args>(args)...)),
      stepper_( ::pressio::ode::create_implicit_stepper<TotalNumberOfDesiredStates>(*rjPolicyOrFullyDiscreteSystem_))
  {}

  stepper_type & lspgStepper(){ return stepper_; }

  template<class SolverType, class ...ArgsOp>
  void operator()(state_type & reducedState,
		  pressio::ode::StepStartAt<independent_variable_type> sStart,
		  pressio::ode::StepCount sCount,
		  pressio::ode::StepSize<independent_variable_type> sSize,
		  SolverType & solver,
		  ArgsOp && ...argsop)
  {
    stepper_(reducedState, sStart, sCount, sSize,
	     solver, std::forward<ArgsOp>(argsop)...);
  }

private:
  std::unique_ptr<FomStatesManager<TrialSubspaceType>> fomStatesManager_;
  std::unique_ptr<ResJacPolicyOrFullyDiscreteSystemType> rjPolicyOrFullyDiscreteSystem_;
  stepper_type stepper_;
};

template<std::size_t TotalNumberOfDesiredStates, class ...Args>
using LspgUnsteadyProblemFullyDiscreteAPI =
  impl::LspgUnsteadyProblem<TotalNumberOfDesiredStates, Args...>;


}}} // end pressio::rom::impl
#endif  // PRESSIO_ROM_IMPL_LSPG_UNSTEADY_PROBLEM_HPP_
