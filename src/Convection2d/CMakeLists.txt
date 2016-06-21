

# add source file
set(SRC_LIST
        ConvectionDriver2d.c
        ConvectionDriver2d.h
        Mesh2d.c
        InitCondition.c
        NcOutput.c
        NcOutput.h
        ConvectionRun2d.c
        ConvectionRHS2d.c
        Postprocess.c)
#        PairFace.c
#        SetUp.c
#        BuildMaps.c
#        InitMeshInfo.c
#        LoadBalance2d.c
#        DisDetector.c
#        LimiterBJ2d.c)

# set executable
add_executable(Convection2d ${SRC_LIST})

# set static library
target_link_libraries(Convection2d pnetcdf)
target_link_libraries(Convection2d parmetis)
target_link_libraries(Convection2d metis)
target_link_libraries(Convection2d PhysDomain)

# install command
install(TARGETS Convection2d
        RUNTIME DESTINATION bin)