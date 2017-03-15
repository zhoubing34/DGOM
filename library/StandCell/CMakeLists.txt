include_directories(${CMAKE_SOURCE_DIR}/library)

set(STDREGIONS_SOURCE
        dg_cell.c
        dg_cell_face.c
        dg_cell_volume.c
        dg_cell_point.c
        dg_cell_line.c
        dg_cell_triangle.c
        dg_cell_quadrilateral.c)

add_library(standcell STATIC ${STDREGIONS_SOURCE})
add_dependencies(standcell utility polylib)
target_link_libraries(standcell utility polylib)

install(TARGETS standcell ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)