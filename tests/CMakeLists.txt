cmake_minimum_required(VERSION 2.6)
# Locate GTest
find_package(GTest)
if (GTEST_FOUND)
    include_directories(${GTEST_INCLUDE_DIRS} ${flat_SOURCE_DIR}/src/)

    add_executable(vec2 vec2.cc)
    target_compile_features(vec2 PRIVATE cxx_aggregate_default_initializers)
    target_link_libraries(vec2 ${GTEST_LIBRARIES} pthread)

    add_executable(rng rng.cc)
    target_compile_features(rng PRIVATE cxx_auto_type)
    target_link_libraries(rng ${GTEST_LIBRARIES} pthread)
endif()
