CMake Keywords
##############

.. important::

   All CMake keywords are prefixed with ``PRESSIO_`` which is case-sensitive.

   Recall that to set a keyword in CMake you used the syntax ``-Dkeyword_name``.


Third-party Libraries (TPLs)
============================

The following options control enabling TPLs:

.. list-table::
   :widths: 30 60 10
   :header-rows: 1
   :align: left

   * - Variable
     - Description
     - Default

   * - ``PRESSIO_ENABLE_TPL_EIGEN``
     - self-explanatory
     - ``ON``

   * - ``PRESSIO_ENABLE_TPL_TRILINOS``
     - self-explanatory
     - ``OFF``

   * - ``PRESSIO_ENABLE_TPL_MPI``
     - self-explanatory
     - ``OFF``  automatically ``ON`` if ``PRESSIO_ENABLE_TPL_TRILINOS=ON``

   * - ``PRESSIO_ENABLE_TPL_KOKKOS``
     - self-explanatory
     - ``OFF``\ ; automatically ``ON`` if ``PRESSIO_ENABLE_TPL_TRILINOS=ON``

   * - ``PRESSIO_ENABLE_TEUCHOS_TIMERS``
     - self-explanatory
     - ``OFF``  automatically ``ON`` if ``PRESSIO_ENABLE_TPL_TRILINOS=ON``


Obviously, the choice of which TPLs to enable is related to
your application's dependency requirements.
For example, if you have an application that relies on
Trilinos data structures and want to use ``pressio-rom``\ ,
then it makes sense to enable the Trilinos dependency.
On the contrary, if you have an application that relies only on
Eigen data structures, then it makes sense to only leave only Eigen on
and disable the rest.

Also, we note that some of the cmake variables listed above are connected
and cannot be turned on individualy.
For example, if we enable Trilinos then ``pressio-rom`` automatically
enables also Kokkos, BLAS, LAPACK and MPI.


Logging
=======

.. list-table::
   :widths: 30 60 10
   :header-rows: 1
   :align: left

   * - Variable
     - Description
     - Default

   * - ``PRESSIO_ENABLE_LOGGING``
     - enables logging with ``pressio-log``
     - ``ON``

   * - ``PRESSIO_ENABLE_COLORIZED_OUTPUT``
     - logging output will be colorized (warnings=yellow, errors=red)
     - ``OFF``

   * - ``PRESSIO_SILENCE_WARNINGS``
     - the logger will omit any warning messages
     - ``OFF``


Other Options
=============

.. list-table::
   :widths: 30 60 10
   :header-rows: 1
   :align: left

   * - Variable
     - Description
     - Default

   * - ``PRESSIO_ENABLE_CXX20``
     - enables C++20 standard
     - ``OFF``; turned on if ``CMAKE_CXX_STANDARD`` is set to 20
