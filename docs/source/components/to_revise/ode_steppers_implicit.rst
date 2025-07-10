.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst

Implicit steppers
=================

Header: ``<pressio/ode_steppers_implicit.hpp>``

Public namespace: ``pressio::ode``


Scope
-----

An "implicit stepper" in pressio is an abstraction that represents
"how" to take a step when applying an :orange:`implicit scheme`
to initial value problems expressable as

.. math::
  :label: ode_implicit_system1

    \frac{d \boldsymbol{y}}{dt} =
    \boldsymbol{f}(\boldsymbol{y},t; ...),  \qquad y(t_0) = y_0

or with a mass matrix:

.. math::
   :label: ode_implicit_system2

    M(\boldsymbol{y}, t, ...) \frac{d \boldsymbol{y}}{dt} =
    \boldsymbol{f}(\boldsymbol{y},t; ...),  \qquad y(t_0) = y_0

or directly in *fully discrete* form as:

.. math::
  :label: ode_implicit_system3

    \boldsymbol{R}(\boldsymbol{y_{n}}, \boldsymbol{y_{n-1}}, ..., t_n, dt_n; ...) = \boldsymbol{0}, \qquad y(t_0) = y_0

where :math:`y` is the state, :math:`f` is the
right hand side (RHS), :math:`t` is the independent variable,
:math:`M` is the mass matrix, and :math:`R` is the residual.
Note that both :math:`f` and :math:`M`
potentially depend on the state and :math:`t`.

.. admonition:: :medium:`Recall the definition of implicit methods`

    :medium:`Implicit methods update the state by solving a (potentially nonlinear) system of equations involving the current, the predicted state and possibly previous states.`


Usage
-----

- if your problem is of the form :eq:`ode_implicit_system1` or  :eq:`ode_implicit_system2`, you can use the following functionalities

  .. toctree::
     :maxdepth: 1

     ode_steppers_implicit_standard_use


- if your problem is of the form :eq:`ode_implicit_system3`, then you can use this:

  .. toctree::
     :maxdepth: 1

     ode_steppers_implicit_fullydiscrete
