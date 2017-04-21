INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library
        ${CMAKE_SOURCE_DIR}/StdRegions)

set(STDREGIONS_SCR
        dg_cell_test_main.c
        point/dg_cell_point_test.c
        line/dg_cell_line_test.c
        triangle/dg_cell_tri_test.c
        quadrilateral/dg_cell_quad_test.c)

add_executable(stand_cell_test ${STDREGIONS_SCR})

add_dependencies(stand_cell_test utility stand_cell)
target_link_libraries(stand_cell_test utility stand_cell)

install(TARGETS stand_cell_test
        RUNTIME DESTINATION test)
