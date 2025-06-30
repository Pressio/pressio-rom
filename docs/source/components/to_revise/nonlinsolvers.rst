.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst

``nonlinear_solvers``
=====================

Header: ``<pressio/solvers_nonlinear.hpp>``

..
   Scope
   -----

   The pressio nonlinear solvers are meant to cover a (currrently) basic
   set of algorithms with the flexibility of adding more if needed, and support arbitrary data types.
   If you are interested in contributing or you need
   other methods, open an issue and/or PR.

   At a high-level, how does this work?
   ------------------------------------

   You are responsible of:

   * A: defining your problem by creating a class with a specific API
   * B: instantiating a solver that is suitable to solve your problem
   * C: (optionally) setting the convergence criteria, tolerances and the update method
   * D: executing the solve

Content
-------

.. toctree::
   :maxdepth: 1

   nonlinsolvers_newton
   nonlinsolvers_gn_neq
   nonlinsolvers_gn_qr
   nonlinsolvers_lm
   nonlinsolvers_enums
   nonlinsolvers_concepts


..
   Currently, we support:

   .. list-table::
      :widths: 30 15 55
      :header-rows: 1

      * - Name
	- Doc Page
	- Usecase
      * - Newton-Raphson
	- `Link <nonlinsolvers_nr.html>`__
	- Systems of nonlinear equations (see `link <https://link.springer.com/content/pdf/bbm%3A978-3-319-69407-8%2F1.pdf>`__ , `link <https://www.cmu.edu/math/undergrad/suami/pdfs/2014_newton_method.pdf>`__ )
      * - Gauss-Newton
	- `Link <nonlinsolvers_gn.html>`__
	- Nonlinear least-squares problem (see `link <https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm>`__ )
      * - Levenberg–Marquardt
	- `Link <nonlinsolvers_lm.html>`__
	- Nonlinear least-squares problem (see `link <https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm>`__ )
      * - Iteratively reweighted least squares (irls)
	- `Link <nonlinsolvers_irls.html>`__
	- Optimization formulated in a p-norm (see `link <https://en.wikipedia.org/wiki/Iteratively_reweighted_least_squares>`__ )

   .. A glimpse of the API
   .. --------------------

   .. To create a solver, we expose specific factory functions for each algorithm.

   .. .. code-block:: cpp

   ..    namespace pnonls = pressio::nonlinearsolvers;
   ..    auto solver = pnonls::create_newton_raphson(....);
   ..    auto solver = pnonls::create_gauss_newton(....);

   .. Please refer to each method's documentation for the details on the API and arguments.


   Practically, how does this work?
   --------------------------------

   You are responsible of:

   * A: defining your problem by creating a class with a specific API
   * B: instantiating a solver that is suitable to solve your problem
   * C: (optionally) setting the convergence criteria, tolerances and the update method
   * D: executing the solve

   Below we provide a brief description of these steps to give you
   an idea of how things work. More details can be found on the individual
   doc page of each method.

   A: :under:`Define the problem`
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   To *define your problem* you need to provide to pressio an
   instance of a "problem class" that you write such that
   it fully encapsulates the computation of the needed operators to operate on.
   Obivously, this problem class must meet a specific API.

   We define four main APIs, namely the "residual-jacobian",
   "fused residual-jacobian", "hessian-gradient" ,
   and "fused hessian-gradient". For the sake of explanation,
   here is a *syntactically* valid residual-jacobian API of a problem:

   .. code-block:: cpp

      struct ProblemResJac
      {
	using state_type     = /* your type */;
	using residual_type  = /* your type */;
	using jacobian_type  = /* your type */;

	state_type    createState() const;
	residual_type createResidual() const;
	jacobian_type createJacobian() const;
	void residual(const state_type& x, residual_type & res) const;
	void jacobian(const state_type& x, jacobian_type & jac) const;
      };


   The key aspect of the class above is that, when instantiated (which obviously
   depends on what you do), it is self-contained and it *fully defines* a problem.
   This is a critical aspect of our API design.
   We want users to specify problems via classes that can completely represent them.
   The above snippet is just syntactic for the sake of explanation.
   More details on the actual concept (syntax + semantics)
   are discussed `here <nonlinsolvers_system_api.html>`_.
   Depending on which concepts your problem meets, you can use
   certain algorithms but not necessarily others.
   The following table helps clarifying the problem API/algorithms admissibility:

   .. list-table::
      :header-rows: 1

      * - Algorithm
	- Residual-Jacobian API
	- Fused Residual-Jacobian API
	- Hessian-Gradient API
	- Fused Hessian-Gradient API
      * - Newton-Raphson
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#10003`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#10003`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -`
      * - Gauss-Newton via Normal-Equations
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#10003`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#10003`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#10003`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#10003`
      * - Gauss-Newton via QR
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#10003`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#10003`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -`
      * - Levenberg-Marquardt
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#10003`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#10003`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#10003`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#10003`
      * - irls
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#10003`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#10003`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -`
	- :raw-html-m2r:`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -`

   Please refer to each method's documentation for the details.


   B: :under:`Instantiate the solver`
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   Once you have the problem, you create a solver using the proper
   the factory functions for the type of solver you want, e.g.:

   .. code-block:: cpp

      namespace pnonls = pressio::nonlinearsolvers;
      auto solver = pnonls::create_newton_raphson(....);
      // ...

   Please refer to each method's documentation for the details on the API and arguments.


   C: :under:`Convergence, tolerance and update`
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   The convergence criterion and associated tolerance are used to decide
   why and when the solver needs to terminate.
   We currently support these termination criteria:

   .. list-table::
      :header-rows: 1

      * - Enum value
	- Applicable to:
      * - ``Stop::AfterMaxIters``
	- all algorithms
      * - ``Stop::WhenCorrectionAbsoluteNormBelowTolerance``
	- all algorithms
      * - ``Stop::WhenCorrectionRelativeNormBelowTolerance``
	- all algorithms
      * - ``Stop::WhenResidualAbsoluteNormBelowTolerance``
	- all algorithms
      * - ``Stop::WhenResidualRelativeNormBelowTolerance``
	- all algorithms
      * - ``Stop::WhenGradientAbsoluteNormBelowTolerance``
	- least-squares solvers
      * - ``Stop::WhenGradientRelativeNormBelowTolerance``
	- least-squares solvers

   which you set/query using the following methods:

   .. code-block:: cpp

      class Solver{
	...

	void setUpdatingCriterion(Update value);
	Update updatingCriterion() const;

	// set/query stopping criterion
	void setStoppingCriterion(Stop value);
	Stop stoppingCriterion() const;

	// set/query max number of iterations
	void setMaxIterations(iteration_type maxIters);
	iteration_type maxIterations() const;

	// this is used to set a single tol for all
	template<class T>
	void setTolerance(T toleranceIn);

	// finer-grained methods for setting tolerances
	void setCorrectionAbsoluteTolerance(/*norm type*/ value);
	void setCorrectionRelativeTolerance(/*norm type*/ value);
	void setResidualAbsoluteTolerance(/*norm type*/ value);
	void setResidualRelativeTolerance(/*norm type*/ value);
	void setGradientAbsoluteTolerance(/*norm type*/ value);
	void setGradientRelativeTolerance(/*norm type*/ value);

	// querying tolerances
	auto correctionAbsoluteTolerance()const;
	auto correctionRelativeTolerance()const;
	auto residualAbsoluteTolerance()const;
	auto residualRelativeTolerance()const;
	auto gradientAbsoluteTolerance()const;
	auto gradientRelativeTolerance()const;
	...
      };

   The update stage represents the *how* the current correction term is combined
   with state to update the latter. We currently support the following:

   .. list-table::
      :widths: 20 25 25 30
      :header-rows: 1

      * - Name
	- Enum value
	- Description
	- Currently supported for:
      * - Default
	- ``Update::Standard``
	- :math:`x_{n+1} = x_{n} + \lambda_{n}`
	- all algorithms
      * - Armijo
	- ``Update::Armijo``
	- :red:`ADD LINK`
	- Gauss-Newton
      * - LM-schedule1
	- ``Update::LMSchedule1``
	- :red:`ADD LINK`
	- Levenberg–Marquardt
      * - LM-schedule2
	- ``Update::LMSchedule2``
	- :red:`ADD LINK`
	- Levenberg–Marquardt

   where :math:`\lambda_{n}` is the correction computed at the n-th iteration of the solver.


   D: :under:`Execute the solve`
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   :red:`SAY SOMETHING`


   |


   Default Settings
   ----------------

   By default, a nonlinear solver uses:

   * update: ``Update::Standard``\ ;
   * stopping: ``Stop::WhenCorrectionAbsoluteNormBelowTolerance``\ ;
   * max number of iterations = 100
   * tolerance = 0.000001 (for everything)

   |


   A note on the solvers' design
   -----------------------------

   If you are interested, here we provide a brief note describing the our design idea.
   The design of the nonlinear solvers has been based on recognizing that, at a very high level,
   a nonlinear solver operates by repeatedly updating a given "state" until a certain criterion is met,
   and each "iteration" involves the following stages:

   *
     A: computing/updating the operators

   *
     B: computing the new correction term

   *
     C: assessing convergence

   *
     D: updating the state using the correction

   This view forms the basis of our design approach: when a solver object is instantiated,
   depending on the chosen algorithm, pressio behind the scenes instantiates the proper
   classes needed for each of the stages above,
   and properly composes them to instantiate the desired solver object.
