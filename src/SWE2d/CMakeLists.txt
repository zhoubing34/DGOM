set(SRC_LIST
        SWEDriver2d.c
        SWEDriver2d.h
        SWEMesh2d.c
        SWEInit2d.c
        SWEOutput2d.c)

add_executable(SWE2d ${SRC_LIST})

target_link_libraries(SWE2d pnetcdf)
target_link_libraries(SWE2d parmetis)
target_link_libraries(SWE2d metis)
target_link_libraries(SWE2d PhysDomain)

# install command
install(TARGETS SWE2d
        RUNTIME DESTINATION bin)