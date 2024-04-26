project(algebrex_test)

file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_executable(${PROJECT_NAME} ${TEST_SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE ${GTEST_INCLUDE_DIR} ${PROJECT_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE GTest::gtest_main)

enable_testing()
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
