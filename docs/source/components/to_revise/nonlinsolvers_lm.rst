
Levenberg-Marquardt
===================

Header: ``<pressio/solvers_nonlinear_levmarq.hpp>``

API
---

.. literalinclude:: ../../../include/pressio/solvers_nonlinear/solvers_create_levenberg_marquardt.hpp
   :language: cpp
   :lines: 60-62, 75-96, 116, 116-117

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
