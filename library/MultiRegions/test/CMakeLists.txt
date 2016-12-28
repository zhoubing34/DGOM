include_directories(${CMAKE_SOURCE_DIR}/library)

set(MULTIREGIONS_TEST_SRC
        mr_test_main.c
        mr_grid_test.c
        mr_grid_test.h
        mr_reg_test.c
        mr_reg_test.h
        mr_mesh_test.c
        mr_mesh_test.h)

add_executable(mulitregions_test ${MULTIREGIONS_TEST_SRC})

add_dependencies(mulitregions_test multiregions)
target_link_libraries(mulitregions_test multiregions)

install(TARGETS mulitregions_test
        RUNTIME DESTINATION test)