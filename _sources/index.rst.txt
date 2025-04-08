.. role:: raw-html-m2r(raw)
   :format: html

.. include:: mydefs.rst

pressio-rom C++ library
===================

Documentation of the `C++ library <https://github.com/Pressio/pressio-rom>`__, one element of the `Pressio ecosystem <https://pressio.github.io/>`_.

This work was started with a focus on projection-based reduced-order models (ROMs),
which is a strongly **multidisciplinary** topic.
Working towards a production-level ROM capability inevitably means spanning
multiple fields ranging from, e.g., linear algebra, nonlinear solvers
and optimization, to time integration, distributed computing and HPC.
This constitutes a substantial challenge to tackle, whose complexity
increases if aiming to develop a **generic** library.

To start such a project from the ground up, grow it and then being
able to maintain it, we believe **modularity, abstractions
and well-defined APIs** to be fundamental design principles to rely on.
This has been, and still is, at the core of our development effort,
and has lead to a highly *modular* design of pressio (see table below):
each component (level) of the stack covers a specific capability and depends,
via *well-defined public APIs*, on the ones below it. This has required (and still does)
a considerable development effort, since each component needs "attention"
and can easily be scoped into an independent, full-time project.

So why doing all this rather than adopting a different, simpler approach, for example,
limiting and hiding as implementation some of the supporting functionalities?
Because we believe the current structure/design offers several major benefits
that would be hard---and in some cases impossible---to obtain otherwise: **flexibility,
extensibility, maintainability, and usability of each component on its own.**
One drawback is that at any point in time, the various components might
have different maturity levels, so reaching a comparable and solid maturity
across the stack might take some time---our current goal is to obtain
in version ``1.0.0`` a uniform maturity level across
the ``rom``, ``ode``, and ``solvers`` components. Please keep this in mind while browsing
the documentation and the code.

The supported capabilities of ``pressio-rom`` are included in the table below.
Each component (level) of the stack depends on the ones below it.

|

.. list-table::
   :widths: 10 48 42
   :header-rows: 1
   :align: left

   * -
     - Description
     - Header(s)

   * - ``rom``
     - concepts :raw-html-m2r:`<br/>` (linear) subspaces :raw-html-m2r:`<br/>` Galerkin: steady :raw-html-m2r:`<br/>` Galerkin: unsteady :raw-html-m2r:`<br/>` LSPG: steady :raw-html-m2r:`<br/>` LSPG: unsteady :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>`
     - ``<pressio/rom_concepts.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_subspaces.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_galerkin_steady.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_galerkin_unsteady.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_lspg_steady.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom_lspg_unsteady.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/rom.hpp>`` :small:`includes all`

   * - ``ode``
     - concepts :raw-html-m2r:`<br/>` explicit steppers :raw-html-m2r:`<br/>` implicit steppers :raw-html-m2r:`<br/>` ``advance_<*>`` fncs :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>`
     - ``<pressio/ode_concepts.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/ode_steppers_explicit.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/ode_steppers_implicit.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/ode_advancers.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/ode.hpp>`` :small:`includes all`

   * - ``solvers_nonlinear``
     - concepts :raw-html-m2r:`<br/>` Newton method :raw-html-m2r:`<br/>` Gauss-Newton :raw-html-m2r:`<br/>` Lev.-Marq. :raw-html-m2r:`<br/>` :raw-html-m2r:`<br/>`
     - ``<pressio/solvers_nonlinear_concepts.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/solvers_nonlinear_newton.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/solvers_nonlinear_gaussnewton.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/solvers_nonlinear_levmarq.hpp>`` :raw-html-m2r:`<br/>` ``<pressio/solvers_nonlinear.hpp>`` :small:`includes all`

   * - ``solvers_linear``
     - linear dense (on-node) solvers
     - ``<pressio/solvers_linear.hpp>``


These components depend on the `Pressio/pressio-ops <https://github.com/Pressio/pressio-ops>`_ header-only library.

Get Started
-----------

* `Configure pressio-rom <configuration.html>`_: (currently) it is a header-only library, so should be trivial

* Explore our `end-to-end ROM demos <https://pressio.github.io/pressio-tutorials/endtoend/readthisfirst.html>`_ to
  see ``pressio-rom`` in action and to experiment directly

* Explore the `tutorials showing the individual capabilities <https://pressio.github.io/pressio-tutorials>`_


Generic programming and concepts
--------------------------------

Arguably the main foundation of pressio-rom is the use of
generic programming--*or, more humbly, we can at least say that it is what we strive for*.
Since the early development stages, we have relied on concept-driven design.
Note, that, if you have used or use C++ templates, you *have* used
concepts knowingly or not. This is because when you write a function or class
template, you have some expectations of what a template needs to expose/do.
C++20 concepts are, in some sense, a way to *explicitly* formalize those expectations.


Until we can stably upgrade to C++20, we cannot by default use C++20 concepts,
so we currently guard the concepts in pressio-rom inside a
preprocessor directive ``#ifdef PRESSIO_ENABLE_CXX20``.
The behavior is as follows:

- if ``PRESSIO_ENABLE_CXX20`` is *enabled*: concepts are compiled and
  enforced *stricto sensu* on the pressio APIs as discussed by this documentation

- if ``PRESSIO_ENABLE_CXX20`` is *disabled*: this is the default case because the
  default pressio C++ standard is currently C++17. In this case, the "C++20 concepts"
  are not compiled but the constraints they represent are still valid and implemented
  differently such that their enforcement is done via a combination of SFINAE and static asserts.

.. important::

   Well-defined concepts are hard to design and it takes time. Concepts used in pressio are
   still being developed. Some are more mature than others. The approach we adopt is to first
   focus on the syntax, then then we will revise them with proper semantics. Keep this in mind
   if some concepts seem incomplete.

Read more about `concepts in Pressio <concepts_note.html>`_.

..
   Here, the term concept does not necessarily
   refer to the C++ concepts feature introduced in C++20.
   You can think of it more broadly as "what properties/syntax a type meets,
   what you can do with it and, also, what a type should definitely satisfy".
   The message we want to convey is that *"concepts" are a fundamental
   design part of pressio*. In our documentation, we make the effort to
   highlight the use of concepts
   by dedicating to each component of the library a full section
   to discuss and formalize how concepts are used in that component.


License and Citation
--------------------

The full license (BSD-3) is available `here <https://github.com/Pressio/pressio/blob/main/LICENSE>`_.

Sooner or later we will publish this... in the meantime, you can find on arXiv
an (outdated) preprint at: https://arxiv.org/abs/2003.07798

Questions?
----------

Find us on Slack: https://pressioteam.slack.com or
open an issue on `github <https://github.com/Pressio/pressio>`_.


.. toctree::
   :maxdepth: 1
   :hidden:

   configuration
   keywords
   concepts_note

.. toctree::
   :caption: API
   :maxdepth: 1
   :hidden:

   ./components/rom
   ./components/ode
   ./components/nonlinsolvers
   ./components/linsolvers

.. toctree::
   :caption: Miscellanea
   :maxdepth: 1
   :hidden:

   Tutorials <https://pressio.github.io/pressio-tutorials>
   GitHub Repo <https://github.com/Pressio/pressio-rom>
   Open an issue/feature req. <https://github.com/Pressio/pressio-rom/issues>
   license
