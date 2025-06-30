
.. include:: ../mydefs.rst

LSPG: Steady
============

Header: ``<pressio/rom_lspg_steady.hpp>``

API
---

.. literalinclude:: ../../../include/pressio/rom/lspg_steady.hpp
   :language: cpp
   :lines: 14-23, 9, 37-48, 9, 64-75, 9, 88-101


..
   Description
   ~~~~~~~~~~~

   Overload set to instantiate a default or hyper-reduced problem problem (1), or masked problem (2).

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


   - `rom::lspg::steady::ComposableIntoDefaultOrHyperReducedProblem <rom_concepts_steady_lspg/default.html>`__

   - `rom::lspg::steady::ComposableIntoMaskedProblem <rom_concepts_steady_lspg/masked.html>`__


   Preconditions
   ~~~~~~~~~~~~~

   .. _steadyLspgPreconditions:

   1. all arguments passed to ``create_steady_problem`` must have a
      lifetime *longer* that that of the instantiated problem, i.e., they must be
      destructed *after* the problem instantiated goes out of scope

   2. the trial subspace must be an admissible approximation
      of the specific full state/problem represented by the ``fomSystem`` instance

   Return value, Postconditions and Side Effects
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   - the return value is an instance of class representing a steady LSPG problem.

       The return type is implementation defined, but guaranteed to
       model the ``SystemWithFusedResidualAndJacobian``
       concept discussed `here <nonlinearsolvers_concepts/c2.html>`__.

   - Purely syntactically, the problem class API is:

     .. code-block:: cpp

	 // This is not the actual class, it just describes the API
	 template<class TrialSubspaceType, ...> // exposition only
	 class SteadyLspgProblemExpositionOnly
	 {
	   public:
	     using state_type     = typename TrialSubspaceType::reduced_state_type;
	     using residual_type  = /*see description below*/;
	     using jacobian_type  = /*see description below*/;

	     state_type    createState() const;
	     residual_type createResidual() const;
	     jacobian_type createJacobian() const;
	     void residualAndJacobian(const state_type & reducedState,
				      residual_type & R,
				      jacobian_type & J,
				      bool computeJacobian) const;
	 };

     where:

     - ``state_type`` aliases the reduced state type of your trial subspace class

     - ``residual_type``, ``jacobian_type`` are as follows:

       - ``std::is_same<residual_type, typename FomSystemType::residual_type>::value == true``

       - ``std::is_same<jacobian_type,
	 decltype(std::declval<FomSystemType const &>().createResultOfJacobianActionOn(trialSubspace.basis()))>::value == true``

   - any necessary memory allocation needed for the implementation
     occurs when the constructor of the class is called. However, we
     guarantee (for now) that the implementation only uses via *const references*
     (as opposed to copying) the arguments passed to the ``create_steady_problem``.
     This is why it is critical to ensure :ref:`precondition 1 <steadyGalerkinPreconditions>`
     is satisfied.

   Solve the problem
   -----------------

   Solving a steady LSPG problem boils down to solving
   a (nonlinear) least-squares problem.
   The problem class exposes the operators defining the system to solve.

   Using pressio nonlinear solvers
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   :red:`need to say more`

   .. important::

      The code below only shows what is needed when using data types
      for which pressio already has specializations inside pressio/ops.
      In this case, the pressio nonlinear solvers will compile and work
      without any additional information.
      If you use data types that are "custom", you also need to specialize
      any ops needed by the solver you choose, so the following code won't compile.

   .. code-block:: cpp

      int main()
      {
	namespace pls   = pressio::linearsolvers;
	namespace pnls  = pressio::nonlinearsolvers;
	namespace plspg = pressio::rom::lspg;

	// assuming trialSubspace and fomSystem already exist
	auto problem = pgal::create_steady_problem(trialSubspace, fomSystem);

	// let's say we want to use Gauss-newton, then we can do this
	using reduced_state_type = typename decltype(problem)::state_type;
	using default_types = pressio::rom::SteadyLspgDefaultOperatorsTraits<reduced_state_type>;
	using gradient_type = typename default_types::gradient_type;
	using hessian_type  = typename default_types::hessian_type;

	using solver_tag = pls::direct::HouseholderQR;
	using linear_solver_t = pls::Solver<solver_tag, hessian_type>;
	auto nonLinearSolver  = pnls::create_gauss_newton(problem, linear_solver_t{});

	auto reducedState = problem.createState();
	// set initial condition for reducedState somehow
	// set other parameters for the solver if needed
	nonLinearSolver.solve(problem, reducedState);
      }


   Use your own solver
   ~~~~~~~~~~~~~~~~~~~

   :red:`need to say more`

   If you don't want to use the pressio solvers,
   you can set up your own because the problem object
   fully identify the system to solve.

   .. code-block:: cpp

      class CustomSolver{
	// constructor as needed

	template<class T>
	void doSolve(const T & problem,
		     typename T::state_type & state)
	{
	  auto R = problem.createResidual();
	  auto J = problem.createJacobian();

	  // note: R and J represent the operators of a (possibly nonlinear)
	  // least-squares problem so that you have to solve accordingly
	  // To evaluate the operators for the given "state":
	  problem.residualAndJacobian(state, R, J, true);
	}
      };

      int main()
      {
	namespace plspg = pressio::rom::lspg;

	// assuming trialSubspace and fomSystem already exist
	auto problem = pgal::create_steady_problem(trialSubspace, fomSystem);

	CustomSolver mySolver;
	auto reducedState = problem.createState();
	// set initial condition for reducedState somehow
	mySolver.doSolve(problem, reducedState);
      }
