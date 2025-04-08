.. role:: raw-html-m2r(raw)
   :format: html

Configuration and Dependencies
==============================

.. tip::

    pressio-rom is header-only, so it does not need to be precompiled and linked to.

.. warning::

    To use pressio-rom, you need at least a C++17 compiler.

Dependencies
------------

Some parts of ``pressio-rom`` contain code and implementations
that are specific to third-party libraries (TPLs).
The main reason for doing this is that we aim, where possible,
to alleviate the user from writing custom operations and allow ``pressio-rom`` to decide when and how to leverage
the native libraries' operations to obtain the best performance.
This should facilitate the integration and use of ``pressio-rom`` by existing applications.
This is a growing capability and we currently only
provide built-in support to some external HPC libraries (see below).
Obviously, these TPLs-specific specializations need to be guarded with
preprecessor directives, and enabled only if one can access the TPLs.

Enabling/disabling specific dependencies can be done by
defining the macros `listed here <keywords.html>`__.
When building the tests, these macros can also be set during
configuration with cmake.


Optional vs Required
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 10 50 40
   :align: left

   * - TPL Name
     - Optional/Required
     - Version Known to Work/run in CI
   * - Eigen
     - Required
     - 3.4.0
   * - Trilinos
     - Optional
     - commits: 0dc4553, 5bbda25
   * - MPI
     - Optional
     - --
   * - Kokkos
     - Optional
     - 4.4.01
   * - BLAS
     - Optional
     - --
   * - LAPACK
     - Optional
     - --
   * - GoogleTest
     - Optional
     - 1.14.0

Eigen is the only required dependency because it is the
default choice for instantiating the ROM data structures
and solving the (dense) ROM problem.

In practice, what are the steps to get and use Pressio?
-------------------------------------------------------

.. TODO: Add spack instructions once merged

1. Clone `pressio-rom <https://github.com/Pressio/pressio-rom>`_ (defaults to the main branch),
or you can pick a `release version <https://github.com/Pressio/pressio-rom/releases>`_

2. Clone `pressio-ops <https://github.com/Pressio/pressio-ops>`_ and, optionally, `pressio-log <https://github.com/Pressio/pressio-log>`_

3. Add the ``pressio-rom/include``, ``pressio-ops/include``, and (optionally) ``pressio-log/include`` to your project's include directories.

4. Use the Pressio ecosystem

Before including any Pressio files, be sure to identify any third-party dependencies as discussed above.
You can define macros to enable or disable code related to these dependencies.

The following code enables the ``MPI`` and ``Kokkos`` libraries, and activates the Pressio logger.

.. code-block:: cpp

    #define PRESSIO_ENABLE_TPL_MPI 1
    #define PRESSIO_ENABLE_LOGGER 1
    #define PRESSIO_ENABLE_TPL_KOKKOS 1

    #include "pressio/what_you_need.hpp"
    // ...
    int main() {
     // do something
    }

When enabling TPLs with the macros above, pressio
does not need to know where a target TPL exists in your system.
By setting the macro, you are simply telling Pressio that you have
a certain TPL and want to enable the corresponding code in pressio.
The TPLs will be needed at linking stage when you build *your* code that *uses* pressio.

.. tip::

    More information regarding the Pressio logger, including macros and configuration
    details, can be found in the `README <https://github.com/Pressio/pressio-log>`_.

Testing
-------

The following steps explain how to build and runs the Pressio tests.

1. Begin by cloning `pressio-rom <https://github.com/Pressio/pressio-rom>`_.

.. tip::

  You do not need to explicitly clone ``pressio-ops``, as it will be included
  automatically when the tests are built.

  Similarly, ``pressio-log`` is included automatically if the ``PRESSIO_ENABLE_LOGGING``
  CMake variable is turned on.

2. Create a build directory.

.. code-block:: bash

   cd <where-you-cloned-pressio-rom>
   mkdir build && mkdir install

3. Use CMake to configure by passing to the command line the target list of CMake variables to define.

For example, suppose we want to enable support for Trilinos and the logger. We would do:

.. code-block:: bash

   export PRESSIO_ROM_SRC=<where-you-cloned-pressio-rom>
   cd ${PRESSIO_ROM_SRC}/build

   cmake -D PRESSIO_ENABLE_TPL_TRILINOS=ON \
         -D PRESSIO_ENABLE_LOGGING=ON \
         -D PRESSIO_ENABLE_TESTS=ON \
         ${PRESSIO_ROM_SRC}

   make # tests are built

Note that this step **only builds tests** because ``pressio-rom`` is header-only.

By default, this step will also clone and link to the ``Pressio/pressio-ops`` library,
which contains essential operations for ``pressio-ROM``.

.. tip::

  Since the tests assume the role of an application using pressio-rom, they will need
  to link against any TPLs that you enable. To specify the location of a library,
  use the following CMake variable: ``-D <tpl>_DIR=/path/to/tpl/install``.

4. Run the tests

.. code-block:: bash

  cd <where-you-cloned-pressio-rom>/build
  ctest -j <np>
