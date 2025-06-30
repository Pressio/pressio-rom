Newton-Raphson
==============

Header: ``<pressio/solvers_nonlinear_newton.hpp>``

API
---

.. literalinclude:: ../../../include/pressio/solvers_nonlinear/solvers_create_newton.hpp
   :language: cpp
   :lines: 60-62, 80-110, 130-131


Parameters
~~~~~~~~~~

.. list-table::
   :widths: 18 82
   :header-rows: 1
   :align: left

   * -
     -

   * - ``system``
     - instance of your problem

   * - ``linSolver``
     - linear solver to use within each nonlinear iteration

Constraints
~~~~~~~~~~~

Concepts are documented `here <nonlinsolvers_concepts.html>`__.
Note: constraints are enforced via proper C++20 concepts when ``PRESSIO_ENABLE_CXX20`` is enabled,
otherwise via SFINAE and static asserts.

Examples
--------

.. admonition:: Demos
   :class: tip

   1. `full demo <https://pressio.github.io/pressio-tutorials/using_eigen/nonlinsolvers1.html>`__
