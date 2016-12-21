include_directories(${CMAKE_SOURCE_DIR}/library)

set(MULTI_REGIONS_SOURCE
        mr_grid_loadBalance2d.c
        mr_grid_loadBalance2d.h
        mr_grid.c
        mr_grid.h
        mr_grid_resortEToV.c
        mr_grid_resortEToV.h
        mr_grid_uniformGrid.c
        mr_grid_uniformGrid.h
        mr_reg.c
        mr_reg.h
        mr_mesh.c
        mr_mesh.h
        mr_mesh_cellConnect.c
        mr_mesh_cellConnect.h
        mr_mesh_parallelPairs.c
        mr_mesh_parallelPairs.h
        mr_mesh_nodeConnect.c
        mr_mesh_nodeConnect.h)

add_library(MultiRegions STATIC ${MULTI_REGIONS_SOURCE})
add_dependencies(MultiRegions StandCell)
target_link_libraries(MultiRegions
        LibUtilities
        StandCell
        parmetis
        metis)

install(TARGETS MultiRegions ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)