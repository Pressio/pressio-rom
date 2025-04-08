.. role:: raw-html-m2r(raw)
   :format: html

.. include:: mydefs.rst

Concepts in Pressio
===================

Concept-driven design has been and still is critical in
the development of pressio. Here we describe some fundamental
notions that are at the core of our design approach.


``concept = syntactic constraints + semantics``
-----------------------------------------------

In pressio-rom, we follow the approach in "Design of Concept Libraries for C++" (2011)
in defining a concept to be a *meaningful combination of syntactic requirements + axioms/semantics*.
Pure syntax can be fully specified and checked by compilers supporting C++20.
On the contrary, axioms/semantics properties cannot be spelled out and statically checked.
Axioms are, in fact, invariants *assumed* to hold. In the C++ standard, axioms
and semantics are specified/listed as comments. Here, we adopt a different
approach to present axioms (similar to "Design of Concept Libraries for C++", 2011),
that we believe benefits readability: axioms are explicitly listed using
the keyword "axiom". Note, however, that if you want to compile the concepts below
using a C++20 compiler, you have to comment out all axioms.



- :red:`concepts are hard and it will take a few iterations to get meaninngful concepts done`

- :red:`not all components have the same maturity, so this will be reflected in the corresponding concepts too`

- :red:`the documentation of each concept follows the standard: we have syntax plus semantics`


"Satisfying" vs. "modeling" a concept
-------------------------------------

`Quoting the C++ standard <https://en.cppreference.com/w/cpp/concepts>`__: "It is said
that a standard concept is satisfied if its syntactic requirements
are met, and is modeled if it is satisfied and its
semantic requirements (if any) are also met."



..
   The layered design of concepts in pressio
   -----------------------------------------

   If you read the section above (hopefully you have!), you now know that,
   by definition, we interpret concepts are a combination of constraints plus semantics.
   Typically, designing concepts is a multi-step process. One good starting point
   could be to specify the "syntactic" requirements. Then, one adds to this semantics
   requirements and, if needed, revises the syntax, and so on.
   This design process is useful because it allows one to create concepts in a
   gradual/progressive fashion, as if one were adding one layer at a time.
   In pressio, this is exactly what has happened. Some of our deepest concepts
   originally started as just syntactic requirements, and then they were improved
   by adding more layers. In the documentation, you will notice that some concept
   are intentially presented to reflect this layered structure.
   More specifically, you will notice that we have designed our concepts
   such that they reflect the following layered structure:

   .. code:: none

      concepts = constraints + execution/language_axioms + math_axioms

   The execution axioms are meant to cover requirements related to the
   code execution or other language-related details.
   The math axioms instead are meant to cover the mathematical foundations.
   Note that not all concepts have both execution and math axioms.

   .. warning::

      As stated on the main doc page, until we can upgrade to C++20,
      we cannot officially use them *stricto sensu*.
      Concepts are thus currently enforced in the code via SFINAE.
