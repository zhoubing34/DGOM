include_directories(${CMAKE_SOURCE_DIR}/library)

set(UTILITIES_TEST_SOURCE
        utility_test.c
        utility_test.h
        MatrixInverse_test.c
        MatrixMultiply_test.c)

add_executable(utility_test ${UTILITIES_TEST_SOURCE})
add_dependencies(utility_test utility)
target_link_libraries(utility_test utility)

install(TARGETS utility_test
        RUNTIME DESTINATION test)