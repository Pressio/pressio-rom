.. role:: raw-html-m2r(raw)
   :format: html

.. include:: ../mydefs.rst

``advance_n_steps``
===================

Header: ``<pressio/ode_advancers.hpp>``


API
---

.. code-block:: cpp

  namespace pressio { namespace ode{

  template<
    class StepperType,
    class StateType,
    class IndVarType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires ExplicitStepper<StepperType>
  #endif
  void advance_n_steps(StepperType & stepper,                        (1)
		       StateType & state,
		       const IndVarType & startVal,
		       const IndVarType & stepSize,
		       ::pressio::ode::StepCount numSteps);

  template<
    class StepperType,
    class StateType,
    class IndVarType,
    class StepSizePolicyType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires ExplicitStepper<StepperType>
	  && StepSizePolicy<StepSizePolicyType, IndVarType>
  #endif
  void advance_n_steps(StepperType & stepper,                        (2)
		       StateType & state,
		       const IndVarType & startVal,
		       const StepSizePolicyType & stepSizePolicy,
		       ::pressio::ode::StepCount numSteps);

  template<
    class StepperType,
    class StateType,
    class IndVarType,
    class ObserverType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires ExplicitStepper<StepperType>
	  && StateObserver<ObserverType,IndVarType, StateType>
  #endif
  void advance_n_steps(StepperType & stepper,                        (3)
		       StateType & state,
		       const IndVarType & startVal,
		       const IndVarType & stepSize,
		       ::pressio::ode::StepCount numSteps,
		       ObserverType && observer);


    class StepperType,
    class StateType,
    class IndVarType,
    class StepSizePolicyType,
    class ObserverType>
  #ifdef PRESSIO_ENABLE_CXX20
    requires ExplicitStepper<StepperType>
	  && StepSizePolicy<StepSizePolicyType, IndVarType>
	  && StateObserver<ObserverType, IndVarType, StateType>
  #endif
  void advance_n_steps(StepperType & stepper,                        (4)
		       StateType & state,
		       const IndVarType & startVal,
		       const StepSizePolicyType & stepSizePolicy,
		       ::pressio::ode::StepCount numSteps,
		       ObserverType && observer);

  template<
    class StepperType,
    class StateType,
    class IndVarType,
    class AuxT,
    class ...Args>
  #ifdef PRESSIO_ENABLE_CXX20
    requires ImplicitStepper<StepperType, AuxT, Args...>
  #endif
  void advance_n_steps(StepperType & stepper,                        (5)
		       StateType & state,
		       const IndVarType & startVal,
		       const IndVarType & stepSize,
		       ::pressio::ode::StepCount numSteps,
		       AuxT && auxArg,
		       Args && ... args);

  template<
    class StepperType,
    class StateType,
    class StepSizePolicyType,
    class IndVarType,
    class AuxT,
    class ...Args>
  #ifdef PRESSIO_ENABLE_CXX20
    requires ImplicitStepper<StepperType, AuxT, Args...>
	  && StepSizePolicy<StepSizePolicyType, IndVarType>
  #endif
  void advance_n_steps(StepperType & stepper,                        (6)
		       StateType & state,
		       const IndVarType & startVal,
		       StepSizePolicyType && stepSizePolicy,
		       ::pressio::ode::StepCount numSteps,
		       AuxT && auxArg,
		       Args && ... args);

  template<
    class StepperType,
    class StateType,
    class ObserverType,
    class IndVarType,
    class AuxT,
    class ...Args>
  #ifdef PRESSIO_ENABLE_CXX20
    requires ImplicitStepper<StepperType, AuxT, Args...>
	  && StateObserver<ObserverType, IndVarType, StateType>
  #endif
  void advance_n_steps(StepperType & stepper,                        (7)
		       StateType & state,
		       const IndVarType & startVal,
		       const IndVarType & stepSize,
		       ::pressio::ode::StepCount numSteps,
		       ObserverType && observer,
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
    requires ImplicitStepper<StepperType, AuxT, Args...>
	  && StepSizePolicy<StepSizePolicyType, IndVarType>
	  && StateObserver<ObserverType, IndVarType, StateType>
  #endif
  void advance_n_steps(StepperType & stepper,                        (8)
		       StateType & state,
		       const IndVarType & startVal,
		       StepSizePolicyType && stepSizePolicy,
		       ::pressio::ode::StepCount numSteps,
		       ObserverType && observer,
		       AuxT && auxArg,
		       Args && ... args);

  }} // end namespace pressio::ode


Description
-----------

Perform :math:`n` steps of a stepper object and update a state.


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

   * - ``startVal``
     - starting value of the independent variable

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









.. *
..   ``args...``\ : optional objects to forward to stepper's ``operator()`` that are
..   potentially needed to perform one step.
..   Note that these are optional, because your stepper might not need anything.
..   The advance functions will simply forward all these to the ``operator()``
..   of the stepper object. Is having these optional argument really needed?
..   One might argue that if the stepper needs to have access to
..   some objects, then these auxiliary objects can be passed to the stepper
..   upon construction, and so the stepper would already have access to them
..   when its ``operator()`` is called.
..   This is true *if* there is a valid reason for making these auxiliary object
..   data members of the ExplicitStepper class.
..   The problem with this approach is that the stepper would need to know
..   about all theses other types, and this creates a very tight
..   coupling between different objects,
..   that might be not even needed.
..   By letting users provide a pack to the advance functions,
..   we are able to decouple this structure such that any object
..   *only* needed within the stepper's ``operator()`` can be directly passed.
..   To explain this, look at the example below.

.. |

.. When are the variadic arguments useful?
.. ---------------------------------------

.. For example in a scenario like the following:

.. .. code-block:: cpp

..    struct ExplicitStepperClass
..    {
..      template<class AuxiliaryType>
..      void operator()(StateType & state,
..                      const TimeType current_time,
..                      const TimeType time_step_size_to_use,
..                      const int32_t step_count,
..                      AuxiliaryType & aux)
..      {
..         const int value = aux.doComplicateCalculation();
..         if (value % 2 == 0){
..           // update state somehow
..         }
..         else{
..           // update state in a different way
..         }
..      }
..    };

..    class Foo{
..      int doComplicateCalculation(){}
..    };

..    class Bar{
..      int doComplicateCalculation(){}
..    };

..    template<class AuxType>
..    void run()
..    {
..      AuxType a;
..      advance_n_steps(stepper, /*state*/, /*start t*/, /*step size*/, /*num steps*/, a);
..    }

..    int main()
..    {
..      ExplicitStepperClass stepper;
..      run<Foo>();
..      run<Bar>();
..    }

.. Here, we don't want to parametrize the ``StepperClass``
.. on the auxiliary class needed because we don't need it as member.
.. We just need to access the auxiliary class inside the ``operator()``.
