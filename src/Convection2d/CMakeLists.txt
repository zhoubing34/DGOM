# add source file
set(SRC_LIST
        ConvectionDriver2d.c
        ConvectionDriver2d.h
        ConvectionMesh2d.c
        InitCondition.c
        ConvectionOutput2d.c
        ConvectionRun2d.c
        ConvectionRHS2d.c
        Postprocess.c)

# set executable
add_executable(Convection2d ${SRC_LIST})

# set static library
target_link_libraries(Convection2d pnetcdf)
target_link_libraries(Convection2d parmetis)
target_link_libraries(Convection2d metis)
target_link_libraries(Convection2d PhysDomain)

# install command
install(TARGETS Convection2d RUNTIME DESTINATION bin)