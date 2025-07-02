
``advance_*`` functions
=======================

.. toctree::
    :maxdepth: 1

    ode_advance_n_steps
    ode_advance_to_target_point
    ode_advance_to_target_point_with_step_recovery


.. `Free functions implementing various strategies
.. to "advance/update" some information (i.e., stored as a "state")
.. using a "steppable" object (more details below).`

.. Why is this useful? Suppose that you have a generic application/usecase
.. involving the concept of a discrete time or just discrete steps,
.. a notion of "state" that defines your information at a specific step,
.. and a "stepper" object, i.e. an object/abstraction that accepts a state
.. and knows how to take a "step" to update that state.
.. In such a scenario, it can be useful to have different strategies
.. to advance the state, for example: advancing for a fixed number of times,
.. advancing until a certain condition is met, etc.
.. This is what the ``pressio/ode_advancers`` provide.

.. Note this is not specific to applied math or scientific computing applications,
.. but can be something generic. Restricting the context to scientific computing,
.. however, one immediately recognizes these functions to be useful
.. for *time integration* of dynamical systems.
.. This is the reason why these functionalities are currently inside ``pressio::ode``\ ,
.. and this is also why we specifically use "time" below.
.. Later on, we might move or generalize them.

.. The set of capabilities will be extended.
.. If you need something but it is not here, open an issue or
.. a PR on `github <https://github.com/Pressio/pressio>`_.
