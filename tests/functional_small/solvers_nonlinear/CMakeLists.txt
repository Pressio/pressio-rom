set(ROOTNAME ${TESTING_LEVEL}_solvers_nonlinear)

if (PRESSIO_ENABLE_TPL_EIGEN)
  set(name doc_main1)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")

  set(name doc_main2)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")
endif()


# -----------------------------
# newton method
# -----------------------------
if (PRESSIO_ENABLE_TPL_EIGEN)
  set(name newton_problem1_eigen)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc)
  add_serial_utest(${TESTING_LEVEL}_solvers_nonlinear_${name} ${SRC1})

  set(name newton_problem1_eigen_fail)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc)
  add_serial_utest(${TESTING_LEVEL}_solvers_nonlinear_${name} ${SRC1})

  set(name newton_matrixfree_problem1_eigen)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc)
  add_serial_utest(${TESTING_LEVEL}_solvers_nonlinear_${name} ${SRC1})

  set(name newton_problem2_eigen)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc)
  add_serial_utest(${TESTING_LEVEL}_solvers_nonlinear_${name} ${SRC1})

  set(name newton_custom_types_compile_only)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc)
  add_serial_utest(${TESTING_LEVEL}_solvers_nonlinear_${name} ${SRC1})
endif()

# -----------------------------
# gauss-newton with normal equations
# -----------------------------
if (PRESSIO_ENABLE_TPL_EIGEN)
  set(name gaussnewton_normaleqs_problem3_eigen)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")

  set(name gaussnewton_normaleqs_problem5_eigen)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")

  set(name gaussnewton_normaleqs_problem9_eigen)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")

  set(name gaussnewton_normaleqs_custom_types_compile_only)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")
endif()

# -----------------------------
# lev-mar with normal equations
# -----------------------------
if (PRESSIO_ENABLE_TPL_EIGEN)
  set(name levmar_normaleqs_problem9_eigen)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")

  set(name levmar_normaleqs_custom_types_compile_only)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")
endif()

# -----------------------------
# weighted gauss-newton or weighted lev-mar with normal equations
# -----------------------------
if (PRESSIO_ENABLE_TPL_EIGEN)
  set(name weighted_gaussnewton_normaleqs_problem3_eigen)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")

  set(name weighted_gaussnewton_normaleqs_custom_types_compile_only)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")
endif()

# -----------------------------
# compact weighted gauss-newton with normal equations
# -----------------------------
if (PRESSIO_ENABLE_TPL_EIGEN)
  set(name compact_weighted_gaussnewton_normaleqs_problem3_eigen)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")

  set(name compact_weighted_gaussnewton_normaleqs_nontrivial_problem3_eigen)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")

  set(name compact_weighted_gaussnewton_normaleqs_custom_types_compile_only)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")
endif()

# -----------------------------
# gauss-newton via QR
# -----------------------------
if (PRESSIO_ENABLE_TPL_EIGEN)
  set(name gaussnewton_qr_problem3_eigen)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")

  set(name gaussnewton_qr_problem8_eigen)
  add_serial_exe_and_test(${name} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${name}.cc "PASSED")
endif()


# -----------------------------
# add tests for QR
# -----------------------------
add_subdirectory(qr)
