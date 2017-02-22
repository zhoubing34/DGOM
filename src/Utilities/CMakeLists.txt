#
add_executable(dg_gmsh2triangle dg_gmsh2triangle.c)
target_link_libraries(dg_gmsh2triangle utility)
install(TARGETS dg_gmsh2triangle
        RUNTIME DESTINATION bin)

#
add_executable(dg_obc2nc dg_obc2nc.c)
target_link_libraries(dg_obc2nc utility)
install(TARGETS dg_obc2nc
        RUNTIME DESTINATION bin)