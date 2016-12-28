include_directories(${CMAKE_SOURCE_DIR}/library)

set(LIBUTILITIES_TEST_SOURCE
        LibUtilities_test.c
        LibUtilities_test.h
        MatrixInverse_test.c
        MatrixMultiply_test.c)

add_executable(libutility_test ${LIBUTILITIES_TEST_SOURCE})
add_dependencies(libutility_test libutility)
target_link_libraries(libutility_test libutility)

install(TARGETS libutility_test
        RUNTIME DESTINATION test)