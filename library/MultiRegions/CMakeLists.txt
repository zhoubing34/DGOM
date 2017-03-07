include_directories(${CMAKE_SOURCE_DIR}/library)

set(MULTI_REGIONS_SOURCE
        Grid/dg_grid_loadBalance.c
        Grid/dg_grid.c
        Grid/dg_grid_retreatEToV.c
        Grid/dg_grid_reader.c
        Region/mr_reg.c
        Region/dg_reg_volumInfo.c
        Region/dg_reg_surfInfo.c
        Mesh/mr_mesh.c
        Mesh/mr_mesh_cellConnect.c
        Mesh/mr_mesh_parallelPairs.c
        Mesh/mr_mesh_nodeConnect.c
        Mesh/mr_mesh_bc.c)

add_library(multiregions STATIC ${MULTI_REGIONS_SOURCE})
add_dependencies(multiregions standcell)
target_link_libraries(multiregions standcell parmetis metis)

install(TARGETS multiregions ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)