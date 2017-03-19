cmake_minimum_required(VERSION 2.6)

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# Link runTests with what we want to test and the GTest and pthread library
add_executable(runTests src/vec2_test.cc)
target_compile_features(runTests PRIVATE cxx_aggregate_default_initializers)
target_link_libraries(runTests ${GTEST_LIBRARIES} pthread)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -g")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg -g")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Os -DNDEBUG")
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
add_executable(flat src/main.cc src/material.cc)
target_compile_features(flat PRIVATE cxx_aggregate_default_initializers)
