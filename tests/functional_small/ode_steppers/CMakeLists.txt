
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_utest(ode_steppers_concepts ${CMAKE_CURRENT_SOURCE_DIR}/concepts.cc)
endif()


if (PRESSIO_ENABLE_TPL_EIGEN)
  set(name doc_main1)
  add_serial_exe_and_test(
    ${name}
    ode_steppers
    ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")

  set(name doc_main2)
  add_serial_exe_and_test(
    ${name}
    ode_steppers
    ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")

  set(name doc_main3)
  add_serial_exe_and_test(
    ${name}
    ode_steppers
    ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")
endif()


# ========================
#
# EXPLICIT METHODS
#
# ========================
if(PRESSIO_ENABLE_TPL_EIGEN)
  # -----------------
  # NO mass matrix
  # -----------------
  set(FILENAME ode_all_explicit_schemes_no_mass_matrix_rhs_observer_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  # forward euler
  set(FILENAME ode_forward_euler_simple_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  set(FILENAME ode_forward_euler_custom_ind_var_type_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  # rk4
  set(FILENAME ode_rk4_simple_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  set(FILENAME ode_rk4_custom_ind_var_type_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  # ssprk3
  set(FILENAME ode_ssprk3_simple_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  set(FILENAME ode_ssprk3_custom_ind_var_type_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  # ab2
  set(FILENAME ode_ab2_simple_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  set(FILENAME ode_ab2_custom_ind_var_type_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  # -----------------
  # WITH mass matrix
  # -----------------
  set(FILENAME ode_all_explicit_schemes_varying_mass_matrix_rhs_observer_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  # forward euler
  set(FILENAME ode_forward_euler_mass_matrix_manual_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  # with mass matrix, all explicit schemes in sinle exe
  set(FILENAME ode_all_explicit_schemes_varying_mass_matrix_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  set(FILENAME ode_all_explicit_schemes_fixed_mass_matrix_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})
endif()

# ========================
#
# IMPLICIT METHODS (CONT-TIME API)
#
# ========================
if(PRESSIO_ENABLE_TPL_EIGEN)
  # implicit stencil data
  set(FILENAME ode_implicit_stencil_data_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  # bdf1
  set(FILENAME ode_bdf1_simple_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  set(FILENAME ode_bdf1_strong_condition_correctness)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  # with step recovery
  add_serial_exe_and_test(
    integrate_to_time_with_failure_bdf1_stepper
    ode_steppers_implicit
    ode_bdf1_correctness_with_step_recovery.cc "PASSED")

  # bdf2
  set(FILENAME ode_bdf2_simple_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  set(FILENAME ode_bdf2_strong_condition_correctness)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  # cn
  set(FILENAME ode_crank_nicolson_simple_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  # with mass matrix, all implicit schemes in single exe
  set(FILENAME ode_all_implicit_schemes_varying_mass_matrix_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  set(FILENAME ode_all_implicit_schemes_fixed_mass_matrix_correctness_eigen)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  # check calling time is correct
  set(FILENAME ode_all_implicit_schemes_check_app_called_with_correct_time)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})
endif()

# ========================
#
# IMPLICIT METHODS (DISCRETE-TIME API)
#
# ========================
if(PRESSIO_ENABLE_TPL_EIGEN)
  set(FILENAME ode_arbitrary_stepper_check_app_called_with_correct_time)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest(${TESTING_LEVEL}_${FILENAME} ${SRC})

  # arbitrary stepper 1 state, step recovery
  add_serial_exe_and_test(
    integrate_to_time_with_failure_arbitrary_stepper_1_aux_state
    ode_steppers_implicit
    arbitrary_stepper_1_aux_state_with_step_recovery.cc "PASSED")

  add_serial_utest(ode_adjoint_logic_1
    ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_logic_1.cc)
  add_serial_utest(ode_adjoint_logic_2
    ${CMAKE_CURRENT_SOURCE_DIR}/adjoint_logic_2.cc)
endif()


if(PRESSIO_ENABLE_TPL_KOKKOS)
  set(FILENAME ode_implicit_stencil_data_kokkos)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}.cc)
  add_serial_utest_kokkos(${TESTING_LEVEL}_${FILENAME} ${SRC})
endif()
