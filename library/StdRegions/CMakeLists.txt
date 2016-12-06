# =======================================================================================
# StdRegions
# =======================================================================================
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(STDREGIONS_SOURCE
        StdRegions.h
        Triangle.c
        StdRegions.c
        Quadrilateral.c)

add_library(StdRegions STATIC ${STDREGIONS_SOURCE})
add_dependencies(StdRegions LibUtilities Polylib)
target_link_libraries(StdRegions LibUtilities Polylib)


install(TARGETS StdRegions ARCHIVE DESTINATION lib)