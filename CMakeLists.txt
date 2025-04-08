cmake_minimum_required(VERSION 3.18.0 FATAL_ERROR)
project(pressio-rom C CXX)

# tests
#=====================================================================
option(PRESSIO_ENABLE_TESTS "Enable testing for pressio-rom" OFF)
if(PRESSIO_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
