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

add_executable(StandCell_Test ${STDREGIONS_SCR})

add_dependencies(StandCell_Test LibUtilities StandCell)
target_link_libraries(StandCell_Test LibUtilities StandCell)

install(TARGETS StandCell_Test
        RUNTIME DESTINATION test)
