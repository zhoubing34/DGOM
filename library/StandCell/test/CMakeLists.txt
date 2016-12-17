INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library
        ${CMAKE_SOURCE_DIR}/StdRegions)

set(STDREGIONS_SCR
        sc_test.c
        sc_test.h
        sc_tri_test.c
        sc_tri_test.h
        sc_quad_test.c
        sc_quad_test.h)

add_executable(StandCell_Test ${STDREGIONS_SCR})

add_dependencies(StandCell_Test LibUtilities StdRegions)
target_link_libraries(StandCell_Test LibUtilities StdRegions)

install(TARGETS StandCell_Test
        RUNTIME DESTINATION test)
