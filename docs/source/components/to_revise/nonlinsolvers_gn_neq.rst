.. role:: raw-html-m2r(raw)
   :format: html

Gauss-Newton (via normal eqs)
=============================

Header: ``<pressio/solvers_nonlinear_gaussnewton.hpp>``

API
---

.. literalinclude:: ../../../include/pressio/solvers_nonlinear/solvers_create_gauss_newton.hpp
   :language: cpp
   :lines: 60-62, 72-94, 121, 121-155, 233-234

Parameters
~~~~~~~~~~

.. list-table::
   :widths: 18 82
   :header-rows: 1
   :align: left

   * -
     -

   * - ``system``
     - your problem instance

   * - ``linSolver``
     - linear solver to solve the normal equations at each nonlinear iteration

   * - ``weigher``
     - the weighting operator for doing weighted least-squares

Constraints
~~~~~~~~~~~

Concepts are documented `here <nonlinsolvers_concepts.html>`__.
Note: constraints are enforced via proper C++20 concepts when ``PRESSIO_ENABLE_CXX20`` is enabled,
otherwise via SFINAE and static asserts.

Examples
--------

.. admonition:: Demos
   :class: tip

   1. `full demo <https://pressio.github.io/pressio-tutorials/using_eigen/nonlinsolvers2.html>`__


..
   Return Value
   ~~~~~~~~~~~~

   Returns a solver object exposing the following public API:

   .. code-block:: cpp

      // This is not the actual class, it just describes the API
      class Solver
      {

      public:
	// query/set update criterion
	Update currentUpdateCriterion() const;
	void setUpdateCriterion(Update value);

	// query/set stop criterion, tolerance
	Stop currentStopCriterion() const;
	void setStopCriterion(Stop value);
	void setStopTolerance(ScalarType value);
	void setMaxIterations(int newMax);

	void solve(StateType & solutionInOut);
      };
