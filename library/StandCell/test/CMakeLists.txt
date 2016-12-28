INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library
        ${CMAKE_SOURCE_DIR}/StdRegions)

set(STDREGIONS_SCR
        sc_test_main.c
        sc_test_main.h
        triangle/sc_tri_test.c
        triangle/sc_tri_test.h
        triangle/sc_tri_data3.h
        quadrilateral/sc_quad_test.c
        quadrilateral/sc_quad_test.h
        quadrilateral/sc_quad_data3.h)

add_executable(standcell_test ${STDREGIONS_SCR})

add_dependencies(standcell_test libutility standcell)
target_link_libraries(standcell_test libutility standcell)

install(TARGETS standcell_test
        RUNTIME DESTINATION test)
