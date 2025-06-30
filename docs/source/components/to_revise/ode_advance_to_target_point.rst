.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst

``advance_to_target_point``
===========================

Header: ``<pressio/ode_advancers.hpp>``

API
---

.. code-block:: cpp

   template<
     class StepperType,
     class StateType,
     class StepSizePolicyType,
     class IndVarType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires ExplicitStepper<StepperType>
          && StepSizePolicy<StepSizePolicyType, IndVarType>
  #endif
   void advance_to_target_point(StepperType & stepper,          (1)
			       StateType & state,
			       const IndVarType & startVal,
			       const IndVarType & finalVal,
			       StepSizePolicyType && stepSizePolicy);

   template<
     class StepperType,
     class StateType,
     class StepSizePolicyType,
     class ObserverType,
     class IndVarType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires ExplicitStepper<StepperType>
          && StepSizePolicy<StepSizePolicyType, IndVarType>
          && StateObserver<ObserverType, IndVarType, StateType>
  #endif
   void advance_to_target_point(StepperType & stepper,          (2)
			       StateType & state,
			       const IndVarType & startVal,
			       const IndVarType & finalVal,
			       StepSizePolicyType && stepSizePolicy,
			       ObserverType && observer);

   template<
     class StepperType,
     class StateType,
     class StepSizePolicyType,
     class IndVarType,
     class AuxT,
     class ...Args>
  #ifdef PRESSIO_ENABLE_CXX20
    requires ImplicitStepper<StepperType, AuxT, Args...>::value
          && StepSizePolicy<StepSizePolicyType, IndVarType>
          && StateObserver<ObserverType, IndVarType, StateType>
	  && (!StateObserver<AuxT, IndVarType, StateType>)
  #endif
   void advance_to_target_point(StepperType & stepper,          (3)
			       StateType & state,
			       const IndVarType & startVal,
			       const IndVarType & finalVal,
			       StepSizePolicyType && stepSizePolicy,
			       AuxT && auxArg,
			       Args && ... args);

   template<
     class StepperType,
     class StateType,
     class StepSizePolicyType,
     class ObserverType,
     class IndVarType,
     class AuxT,
     class ...Args>
  #ifdef PRESSIO_ENABLE_CXX20
    requires ImplicitStepper<StepperType, AuxT, Args...>::value
          && StepSizePolicy<StepSizePolicyType, IndVarType>
          && StateObserver<ObserverType, IndVarType, StateType>
	  && (!StateObserver<AuxT, IndVarType, StateType>)
  #endif

   void advance_to_target_point(StepperType & stepper,          (4)
			       StateType & state,
			       const IndVarType & startVal,
			       const IndVarType & finalVal,
			       StepSizePolicyType && stepSizePolicy,
			       ObserverType && observer,
			       AuxT && auxArg,
			       Args && ... args);

Description
-----------

Overload set for using a stepper object to update/advance a state
until the independent variable reaches a target value.

Parameters
----------

.. list-table::
   :widths: 18 82
   :header-rows: 1
   :align: left

   * -
     -

   * - ``stepper``
     - object that knows *how to* perform a single step

   * - ``state``
     - the "state" information to update

   * - ``startVal``, ``finalVal``
     - the starting and final value of the independent variable

   * - ``numSteps``
     - how many steps to take

   * - ``stepSizePolicy``
     - functor to set the step size

   * - ``stepSize``
     - *constant* step size to use for each step

   * - ``observer``
     - object to "observe" the state's evolution, which can be used
       to potentially collect necessary data/metrics/statistics or do other things from the state

   * - ``auxArg``, ``args``
     - extra arguments that might be needed by the stepper to execute one step


Constraints
-----------

Concepts are documented `here <ode_concepts.html>`__.
Note: constraints are enforced via proper C++20 concepts when ``PRESSIO_ENABLE_CXX20`` is enabled,
otherwise via SFINAE and static asserts.


Preconditions
-------------

:red:`finish`

Mandates
--------

* ``std::is_same<IndVarType, typename StepperType::independent_variable_type>``

* ``std::is_same<StateType, typename StepperType::state_type>``

Return value
------------

None

Postconditions and Side Effects
-------------------------------

:red:`finish`
