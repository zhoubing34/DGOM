# =======================================================================================
# LocalRegions
# =======================================================================================
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(LOCAL_REGIONS_SOURCE
        LocalRegions.h
        Mapping.c)

add_library(LocalRegions STATIC ${LOCAL_REGIONS_SOURCE})
add_dependencies(LocalRegions LibUtilities Polylib)
target_link_libraries(LocalRegions LibUtilities Polylib)

install(TARGETS LocalRegions ARCHIVE DESTINATION lib)
