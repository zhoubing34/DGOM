include_directories(${CMAKE_SOURCE_DIR}/library)

set(STDREGIONS_SOURCE
        sc_stdcell.h
        sc_stdcell.c
        sc_triangle.c
        sc_quadrilateral.c)

add_library(standcell STATIC ${STDREGIONS_SOURCE})
add_dependencies(standcell libutility polylib)
target_link_libraries(standcell libutility polylib)

install(TARGETS standcell ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)