# This file was generated in part with generative AI

cmake_minimum_required(VERSION 3.16)

set(CATCH_INSTALL_DOCS OFF CACHE BOOL "Do not install Catch2 docs")
set(CATCH_INSTALL_EXTRAS OFF CACHE BOOL "Do not install Catch2 extras")
set(CATCH_BUILD_TESTING OFF CACHE BOOL "Do not build Catch2 self-tests")

if(NOT TARGET Catch2)
add_subdirectory(
  ${PROJECT_SOURCE_DIR}/external/Catch2
  ${CMAKE_BINARY_DIR}/external/Catch2
)
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/external/Catch2/extras)
include(Catch)

set(CHARLES_TEST_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}/src
)

set(CHARLES_UNIT_TEST_SOURCES
  catch2_define.cpp
  test_constants.cpp
  test_sg03.cpp
  test_tqm05.cpp
)

add_executable(charlemagne_unit_tests ${CHARLES_UNIT_TEST_SOURCES})

target_include_directories(charlemagne_unit_tests PRIVATE ${CHARLES_TEST_INCLUDE_DIRS})

target_link_libraries(charlemagne_unit_tests
  PRIVATE
    Catch2::Catch2
    singularity-opac::flags
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_definitions(charlemagne_unit_tests PRIVATE DEBUG=1)
else()
  target_compile_definitions(charlemagne_unit_tests PRIVATE DEBUG=0)
endif()

catch_discover_tests(charlemagne_unit_tests
  DISCOVERY_MODE PRE_TEST
  TEST_PREFIX charlemagne.unit.
)
