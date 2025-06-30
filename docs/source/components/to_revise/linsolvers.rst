``linear_solvers``
==================

Header ``<pressio/solvers_linear.hpp>``

Public namespace: ``pressio::linearsolvers``

Overview
--------

We clearly state upfront that this component of pressio
is not a self-contained, general purpose library of linear solvers
implemented from scratch.
In the future, we might revise this, but for now it is
an **intentional** choice: there exist already many
open-source highly-optimized libraries that provide linear solvers
for both shared-memory, e.g. Armadillo, Eigen, Blaze, Magma, LAPACK,
and distributed environment, e.g. Trilinos, PETSc, HYPRE.


**The primary goal of the pressio's linear solvers is to support usecases that are needed for ROMs.
As discussed in the introduction, ROMs involve "small"
dense systems that suitably fit on a single compute node.
Pressio's linear solvers are thus designed and implemented to be
wrappers (exposing a common interface) to existing shared-memory
linear solvers libraries.**

Currently, we support linear systems using either
Eigen or Kokkos containers. The reason for this is that Eigen
and Kokkos are the supported and preferred choices (for now)
to implement ROMs operators. Kokkos, in particular, allows
to operate both on host and GPUs.
Later on, we can easily extend support to other libraries if needed.
Note that even if the scope of the linear solvers is limited,
we still emphasize that, like all of the others subpackages in pressio,
it can be used idenpendetly.

Why do we need it?
^^^^^^^^^^^^^^^^^^

\todo: explain that we need to a uniform interface.

API
---

.. code-block:: cpp

   pressio::linearsolvers::Solver<Tag, MatrixType>;


* ``Tag``\ : tag type specifying which type of solver to use (more on this below)
* ``MatrixType``\ : data type of the matrix of the linear system to solve
  (used to choose the proper specialization).

We distinguish between iterative and direct solvers.
Obviously, depending on the structure of the matrix, some methods
might not be applicable or might not be efficient.

.. list-table::
   :header-rows: 1

   * - Tag
     - Description
     - Works for:
   * - ``iterative::CG``
     - Conjugate-gradient
     - Eigen
   * - ``iterative::Bicgstab``
     - Biconjugate gradient stabilized
     - Eigen
   * - ``iterative::LSCG``
     - Conjugate-gradient
     - Eigen
   * - ``direct::HouseholderQR``
     - Uses on householder QR
     - Eigen
   * - ``direct::ColPivHouseholderQR``
     - Uses on Householder QR with pivoting
     - Eigen
   * - ``direct::PartialPivLU``
     - Uses LU factorization with partial pivoting
     - Eigen
   * - ``direct::potrsL``
     - Uses Cholesky, lower part
     - Kokkos
   * - ``direct::potrsU``
     - Uses Cholesky, upper part
     - Kokkos
   * - ``direct::getrs``
     - Uses LU factorization
     - Kokkos
   * - ``direct::geqrf``
     - Uses QR fatorization
     - Kokkos

Synopsis
--------

.. code-block:: cpp

   template<class MatrixType>
   class LinearSolver
   {
   public:
     using matrix_type = MatrixType;

     template<class StateType, class RhsType>
     void solve(const MatrixType & A, const RhsType & b, StateType & x);
   };

Example Usage
-------------

\todo: need to discuss difference between iterative and direct

.. code-block:: cpp

   //
   // we want to solve `Ax = b`
   //

   using matrix_type = Eigen::MatrixXd;
   using vector_type = Eigen::VectorXd;
   matrix_type A(10, 10);
   // fill A
   vector_type b(10);
   // fill b

   namespace pls  = pressio::linearsolvers;
   using tag      = pls::direct::HouseholderQr;
   using solver_t = pls::Solver<tag, matrix_type>;
   solver mySolver;
   vector_type x(10);
   mySolver.solve(A, b, x);
