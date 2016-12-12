# =======================================================================================
# MultiRegions
# =======================================================================================
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(MULTI_REGIONS_SOURCE
        MultiRegions.c
        MultiRegions.h
        SetFacePair.c
        SetNodePair.c
        SetElemetPair.c
        LoadBalance.c
        UnstructMesh.c
        UnstructMesh.h
        VertexSort.c
        VertexSort.h)

add_library(MultiRegions STATIC ${MULTI_REGIONS_SOURCE})
add_dependencies(MultiRegions LocalRegions)
target_link_libraries(MultiRegions
        LibUtilities
        LocalRegions
        StdRegions
        parmetis metis)

install(TARGETS MultiRegions ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)