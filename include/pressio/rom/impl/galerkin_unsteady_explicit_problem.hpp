
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_EXPLICIT_PROBLEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_EXPLICIT_PROBLEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

// template <class GalSystem>
// class GalerkinUnsteadyExplicitProblem
// {
//   // note: to deduce the stepper_type it does not really matter
//   // what scheme enum value we use, as long as it is an explicit one
//   using stepper_type =
//     decltype(::pressio::ode::create_explicit_stepper
// 	     (::pressio::ode::StepScheme::ForwardEuler,
// 	      std::declval<GalSystem &>()
// 	      ));

// public:
//   // required aliases to be steppable
//   using state_type = typename GalSystem::state_type;
//   using independent_variable_type  = typename GalSystem::independent_variable_type;

//   template<class ...Args>
//   GalerkinUnsteadyExplicitProblem(::pressio::ode::StepScheme schemeName,
// 				  Args && ... args)
//     : galSystem_(std::forward<Args>(args)...),
//       stepper_( ::pressio::ode::create_explicit_stepper(schemeName, galSystem_) )
//   {}

//   void operator()(state_type & state,
// 		  pressio::ode::StepStartAt<independent_variable_type> sStart,
// 		  pressio::ode::StepCount sCount,
// 		  pressio::ode::StepSize<independent_variable_type> sSize)
//   {
//     stepper_(state, sStart, sCount, sSize);
//   }

// private:
//   GalSystem galSystem_;
//   stepper_type stepper_;
// };

// template <class GalSystem>
// class GalerkinUnsteadyWithMassMatrixExplicitProblem
// {
//   // note: to deduce the stepper_type it does not really matter
//   // what scheme enum value we use, as long as it is an explicit one
//   using stepper_type =
//     decltype(::pressio::ode::create_explicit_stepper
// 	     (::pressio::ode::StepScheme::ForwardEuler,
// 	      std::declval<GalSystem &>()
// 	      ));

// public:
//   // required aliases to be steppable
//   using state_type = typename GalSystem::state_type;
//   using independent_variable_type = typename GalSystem::independent_variable_type;
//   using mass_matrix_type = typename GalSystem::mass_matrix_type;

//   template<class ...Args>
//   GalerkinUnsteadyWithMassMatrixExplicitProblem(::pressio::ode::StepScheme schemeName,
// 						Args && ... args)
//     : galSystem_(std::forward<Args>(args)...),
//       stepper_( ::pressio::ode::create_explicit_stepper(schemeName, galSystem_) )
//   {}

//   template<class LinearSolverType>
//   void operator()(state_type & state,
// 		  pressio::ode::StepStartAt<independent_variable_type> sStart,
// 		  pressio::ode::StepCount sCount,
// 		  pressio::ode::StepSize<independent_variable_type> sSize,
// 		  LinearSolverType & linSolver)
//   {
//     stepper_(state, sStart, sCount, sSize, linSolver);
//   }

// private:
//   GalSystem galSystem_;
//   stepper_type stepper_;
// };

}}} // end pressio::rom::impl
#endif  // PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_EXPLICIT_PROBLEM_HPP_
