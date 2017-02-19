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
        mr_mesh_nodeConnect.h
        mr_mesh_bc.c
        mr_mesh_bc.h)

add_library(multiregions STATIC ${MULTI_REGIONS_SOURCE})
add_dependencies(multiregions standcell)
target_link_libraries(multiregions standcell parmetis metis)

install(TARGETS multiregions ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)