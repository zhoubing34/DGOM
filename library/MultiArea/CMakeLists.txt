include_directories(${CMAKE_SOURCE_DIR}/library)

set(MULTI_REGIONS_SOURCE
        dg_area.c
        Grid/dg_grid_loadBalance.c
        Grid/dg_grid.c
        Grid/dg_grid_retreatEToV.c
        Grid/dg_grid_reader.c
        Grid/dg_grid_connect.c
        Grid/dg_grid_parallelPairs.c
        Grid/dg_grid_BS.c
        Region/dg_region.c
        Region/dg_region_volumInfo.c
        Region/dg_region_surfInfo.c
        Edge/dg_edge.c
        Edge/dg_edge_face_map.c
        Edge/dg_edge_node_map.c
        Edge/dg_edge_surfinfo.c
        Mesh/dg_mesh.c
        Mesh/dg_mesh_connect.c
        Mesh/dg_mesh_fetch_buffer.c
        )

add_library(multi_area STATIC ${MULTI_REGIONS_SOURCE})
add_dependencies(multi_area stand_cell)
target_link_libraries(multi_area stand_cell parmetis metis)

install(TARGETS multi_area ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)