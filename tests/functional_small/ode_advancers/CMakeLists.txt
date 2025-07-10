
set(ROOTNAME ${TESTING_LEVEL}_ode_advancers)

add_serial_utest(
  ${ROOTNAME}_concepts_steppable
  ${CMAKE_CURRENT_SOURCE_DIR}/concepts_steppable.cc)

add_serial_utest(
  ${ROOTNAME}_concepts_step_policy
  ${CMAKE_CURRENT_SOURCE_DIR}/concepts_step_policy.cc)

add_serial_utest(
  ${ROOTNAME}_concepts_state_observer
  ${CMAKE_CURRENT_SOURCE_DIR}/concepts_state_observer.cc)

#
add_serial_utest(
  ${ROOTNAME}_api_check_advance_n_steps
  ${CMAKE_CURRENT_SOURCE_DIR}/api_check_advance_n_steps.cc)

#
add_serial_utest(
  ${ROOTNAME}_advance_n_steps_const_dt
  ${CMAKE_CURRENT_SOURCE_DIR}/advance_n_steps_const_dt.cc)

add_serial_utest(
  ${ROOTNAME}_advance_n_steps_const_dt_with_observer
  ${CMAKE_CURRENT_SOURCE_DIR}/advance_n_steps_const_dt_with_observer.cc)

add_serial_utest(
  ${ROOTNAME}_advance_n_steps_varying_dt
  ${CMAKE_CURRENT_SOURCE_DIR}/advance_n_steps_varying_dt.cc)

add_serial_utest(
  ${ROOTNAME}_advance_n_steps_check_time_and_step
  ${CMAKE_CURRENT_SOURCE_DIR}/advance_n_steps_check_time_and_step.cc)

if(PRESSIO_ENABLE_TPL_EIGEN)
add_serial_exe_and_test(
  advance_to_time_with_failure_mock_stepper
  ${ROOTNAME} to_target_time_with_time_step_recovery.cc "PASSED")
endif()
