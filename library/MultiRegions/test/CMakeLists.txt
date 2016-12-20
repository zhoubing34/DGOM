include_directories(${CMAKE_SOURCE_DIR}/library)

set(MULTIREGIONS_TEST_SRC
        mr_test_main.c
        mr_grid_test.c
        mr_grid_test.h
        mr_reg_test.c
        mr_reg_test.h)

add_executable(MulitRegions_Test ${MULTIREGIONS_TEST_SRC})

add_dependencies(MulitRegions_Test MultiRegions)
target_link_libraries(MulitRegions_Test MultiRegions)

install(TARGETS MulitRegions_Test
        RUNTIME DESTINATION test)