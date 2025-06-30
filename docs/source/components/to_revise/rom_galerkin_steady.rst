
.. include:: ../mydefs.rst

Galerkin: Steady
================

Header: ``<pressio/rom_galerkin_steady.hpp>``

API
---

.. literalinclude:: ../../../include/pressio/rom/galerkin_steady.hpp
   :language: cpp
   :lines: 10-11, 16-25, 40-52, 67-81, 92-93

..
   API
   ---

   .. code-block:: cpp

     namespace pressio { namespace rom{ namespace galerkin{

     template<
       class TrialSubspaceType,
       class FomSystemType>
     #ifdef PRESSIO_ENABLE_CXX20
       requires steady::ComposableIntoDefaultProblem<
		   TrialSubspaceType, FomSystemType>
     #endif
     /*impl defined*/ create_steady_problem(const TrialSubspaceType & trialSubspace,     (1)
					    const FomSystemType & fomSystem);

     template<
       class TrialSubspaceType,
       class FomSystemType,
       class HyperReducerType>
     #ifdef PRESSIO_ENABLE_CXX20
       requires steady::ComposableIntoHyperReducedProblem<
		   TrialSubspaceType, FomSystemType, HyperReducerType>
     #endif
     /*impl defined*/ create_steady_problem(const TrialSubspaceType & trialSubspace,     (2)
					    const FomSystemType & fomSystem,
					    const HyperReducerType & hyperReducer);

     template<
       class TrialSubspaceType,
       class FomSystemType,
       class MaskerType,
       class HyperReducerType>
     #ifdef PRESSIO_ENABLE_CXX20
       requires steady::ComposableIntoHyperReducedMaskedProblem<
		   TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>
     #endif
     /*impl defined*/ create_steady_problem(const TrialSubspaceType & trialSubspace,     (3)
					    const FomSystemType & fomSystem,
					    const maskerType & masker,
					    const HyperReducerType & hyperReducer);

     }}} // end namespace pressio::rom::galerkin

   Description
   ~~~~~~~~~~~

   Overload set to instantiate a default (1), hyper-reduced (2) or masked problem (3).

   Parameters
   ~~~~~~~~~~

   .. list-table::
      :widths: 18 82
      :header-rows: 1
      :align: left

      * -
	-

      * - ``trialSubspace``
	- trial subspace approximating the FOM state space

      * - ``fomSystem``
	- full-order model instance

      * - ``hyperReducer``
	- hyper-reduction operator

      * - ``masker``
	- masking operator

   Constraints
   ~~~~~~~~~~~

   Each overload is associated with a set of constraints.
   With C++20, these would be enforced via concepts using
   the *requires-clause* shown in the API synopsis above.
   Since we cannot yet officially upgrade to C++20, the constraints
   are currently enforced via static asserts (to provide a decent error message)
   and/or SFINAE. The concepts used are:

   - `rom::galerkin::steady::ComposableIntoDefaultProblem <rom_concepts_steady_galerkin/default.html>`__

   - `rom::galerkin::steady::ComposableIntoHyperReducedProblem <rom_concepts_steady_galerkin/hr.html>`__

   - `rom::galerkin::steady::ComposableIntoHyperReducedMaskedProblem <rom_concepts_steady_galerkin/masked.html>`__


   Preconditions
   ~~~~~~~~~~~~~

   .. _steadyGalerkinPreconditions:

   1. all arguments passed to ``create_steady_problem`` must have a
      lifetime *longer* that that of the instantiated problem, i.e., they must be
      destructed *after* the problem instantiated goes out of scope

   2. the trial subspace must be an admissible approximation
      of the specific full state/problem represented by the ``fomSystem`` instance


   Return value, Postconditions, Side Effects
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   - the return is an instance of class representing a Galerkin steady problem

       The return type is implementation defined, but guaranteed
       to model the ``SystemWithFusedResidualAndJacobian``
       concept discussed `here <nonlinearsolvers_concepts/c2.html>`__.


   - Purely syntactically, the problem class API is:

     .. code-block:: cpp

	 // This is not the actual class, it just describes the API
	 template<class TrialSubspaceType, ...> // exposition only
	 class SteadyGalerkinProblemSyntaxOnly
	 {
	   using galerkin_types = SteadyGalerkinTraits<
	     typename TrialSubspaceType::reduced_state_type>;  // exposition only

	   public:
	     using state_type    = typename TrialSubspaceType::reduced_state_type;
	     using residual_type = typename galerkin_types::reduced_residual_type;
	     using jacobian_type = typename galerkin_types::reduced_jacobian_type;

	     state_type    createState() const;
	     residual_type createResidual() const;
	     jacobian_type createJacobian() const;
	     void residualAndJacobian(const state_type & /*reducedState*/,
				      residual_type    & /*reducedResidual*/,
				      jacobian_type    & /*reducedJacobian*/,
				      bool               /*computeJacobian*/) const;
	 };

     where:

     - ``state_type``: must be an Eigen vector as per `this constraints <rom_trial_column_subspace.html>`__

     - ``residual_type``, ``jacobian_type``: types of the *reduced* operators defined as:

       - ``residual_type`` is the same as ``state_type``

       - ``jacobian_type`` is an Eigen dense matrix with column major layout

     - and the following holds:
       ``pressio::all_have_traits_and_same_scalar<state_type, residual_type, jacobian_type>::value == true``


   - any necessary memory allocation needed for the implementation
     occurs when the constructor of the class is called. However, we
     guarantee (for now) that the implementation only uses via *const references*
     (as opposed to copying) the arguments passed to the ``create_steady_problem``.
     This is why it is critical to ensure :ref:`precondition 1 <steadyGalerkinPreconditions>`
     is satisfied.


   Solve the problem
   -----------------

   Solving a steady Galerkin problem practically means solving
   a **reduced dense system of (nonlinear) equations**.
   A Galerkin problem created as shown above exposes the API needed to
   compute the operators defining such reduced system.
   Since the reduced system, by definition, is a *determined system of equations*,
   to solve it one can use a Newton-Raphson solver either
   from the pressio/nonlinear_solvers, or use/implement their own.

   Using pressio nonlinear solvers
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   .. code-block:: cpp

      #include "pressio/rom_subspaces.hpp"
      #include "pressio/rom_galerkin_steady.hpp"

      int main()
      {
	namespace pls  = pressio::linearsolvers;
	namespace pnls = pressio::nonlinearsolvers;
	namespace pgal = pressio::rom::galerkin;

	// assuming trialSubspace and fomSystem already exist
	auto problem = pgal::create_steady_problem(trialSubspace, fomSystem);

	using jacobian_type   = typename decltype(problem)::jacobian_type;
	using linear_solver_t = pls::Solver<pls::iterative::LSCG, jacobian_type>;
	auto nonLinearSolver  = pnls::create_newton_raphson(problem, linear_solver_t{});

	auto reducedState = problem.createState();
	// set initial condition for reducedState somehow
	// set other parameters for the solver if needed
	nonLinearSolver.solve(problem, reducedState);
      }

   Use your own solver
   ~~~~~~~~~~~~~~~~~~~

   If you don't want to use the pressio solvers,
   you can easily set up your own because the problem object
   fully identify the system to solve.

   .. code-block:: cpp

      #include "pressio/rom_subspaces.hpp"
      #include "pressio/rom_galerkin_steady.hpp"

      class CustomSolver{
	// constructor as needed

	template<class T>
	void doSolve(const T & problem,
		     typename T::state_type & state)
	{
	  auto R = problem.createResidual();
	  auto J = problem.createJacobian();

	  for (...){
	    system.residualAndJacobian(state, R, J, true);
	    // do something, update state, etc.
	  }
	}
      };

      int main()
      {
	namespace pgal = pressio::rom::galerkin;

	// assuming trialSubspace and fomSystem already exist
	auto problem = pgal::create_steady_problem(trialSubspace, fomSystem);

	CustomSolver mySolver;
	auto reducedState = problem.createState();
	// set initial condition for reducedState somehow
	mySolver.doSolve(problem, reducedState);
      }
