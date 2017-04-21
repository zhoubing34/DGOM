include_directories(${CMAKE_SOURCE_DIR}/library)

set(STDREGIONS_SOURCE
        dg_cell.c
        dg_cell_face.c
        dg_cell_volume.c
        dg_cell_point.c
        dg_cell_line.c
        dg_cell_triangle.c
        dg_cell_quadrilateral.c)

add_library(stand_cell STATIC ${STDREGIONS_SOURCE})
add_dependencies(stand_cell utility polylib)
target_link_libraries(stand_cell utility polylib)

install(TARGETS stand_cell ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)