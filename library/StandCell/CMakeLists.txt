include_directories(${CMAKE_SOURCE_DIR}/library)

set(STDREGIONS_SOURCE
        sc_stdcell.h
        sc_stdcell.c
        sc_stdcell2d.c
        sc_stdcell2d.h
        sc_triangle.c
        sc_triangle.h
        sc_quadrilateral.c
        sc_quadrilateral.h)

add_library(standcell STATIC ${STDREGIONS_SOURCE})
add_dependencies(standcell utility polylib)
target_link_libraries(standcell utility polylib)

install(TARGETS standcell ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)