cmake_minimum_required(VERSION 2.6)

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Link runTests with what we want to test and the GTest and pthread library
add_executable(runTests vec2_test.cc)
target_compile_features(runTests PRIVATE cxx_aggregate_default_initializers)
target_link_libraries(runTests ${GTEST_LIBRARIES} pthread)

add_executable(scattering main.cc material.cc)
target_compile_features(scattering PRIVATE cxx_aggregate_default_initializers)
