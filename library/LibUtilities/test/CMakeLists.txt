include_directories(${CMAKE_SOURCE_DIR}/library)

set(LIBUTILITIES_TEST_SOURCE
        LibUtilities_test.c
        LibUtilities_test.h
        invM_test.c
        MatrixMultiply.c)

add_executable(LibUtilities_Test ${LIBUTILITIES_TEST_SOURCE})
add_dependencies(LibUtilities_Test LibUtilities)
target_link_libraries(LibUtilities_Test LibUtilities)

install(TARGETS LibUtilities_Test
        RUNTIME DESTINATION test)