set(SRC_LIST
        swe_dirver2d.c
        swe_dirver2d.h
        SWEMesh2d.c
        SWEInit2d.c
        SWEOutput2d.c
        SWERun2d.c
        SWERHS2d.c
        SWEUtilities.c
        swe_getparameter.c
        swe_getparameter.h)

add_executable(swe2d ${SRC_LIST})

target_link_libraries(swe2d pnetcdf)
target_link_libraries(swe2d parmetis)
target_link_libraries(swe2d metis)
target_link_libraries(swe2d physfield)

# install command
install(TARGETS swe2d RUNTIME DESTINATION bin)