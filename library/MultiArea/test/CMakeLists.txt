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

add_executable(multi_area_test ${MULTIREGIONS_TEST_SRC})
add_dependencies(multi_area_test multi_area)
target_link_libraries(multi_area_test multi_area)
install(TARGETS multi_area_test RUNTIME DESTINATION test)