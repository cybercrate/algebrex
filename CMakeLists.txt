cmake_minimum_required(VERSION 3.20)
project(algebrex)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Werror -pedantic")

set(GTEST_ROOT "./lib/googletest")

set(GTEST_INCLUDE_DIR "${GTEST_ROOT}/googletest/include")
set(PROJECT_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

add_subdirectory(${GTEST_ROOT})
add_subdirectory(test)
