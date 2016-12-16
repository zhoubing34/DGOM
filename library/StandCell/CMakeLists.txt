# =======================================================================================
# StdRegions
# =======================================================================================
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(STDREGIONS_SOURCE
        sc_stdcell.h
        sc_stdcell.c
        sc_triangle.c
        sc_quadrilateral.c)

add_library(StdRegions STATIC ${STDREGIONS_SOURCE})
add_dependencies(StdRegions LibUtilities Polylib)
target_link_libraries(StdRegions LibUtilities Polylib)

install(TARGETS StdRegions ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)