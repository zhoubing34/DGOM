include_directories(${CMAKE_SOURCE_DIR}/library)
set(MULTIREGIONS_TEST_SRC
        multiregion_test.c
        multiregion_test.h
        dg_grid_test.c
        dg_grid_test.h
        dg_reg_test.c
        dg_reg_test.h
        dg_mesh_test.c
        dg_mesh_test.h
        dg_edge_test.c)

add_executable(mulitregions_test ${MULTIREGIONS_TEST_SRC})
add_dependencies(mulitregions_test multiregions)
target_link_libraries(mulitregions_test multiregions)
install(TARGETS mulitregions_test
        RUNTIME DESTINATION test)