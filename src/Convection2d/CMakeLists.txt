# add source file
set(SRC_LIST
        ConvectionDriver2d.c Mesh2d.c ../Utilities/Utils.c ../Utilities/ParallelPairs.c
        PairFace.c SetUp.c BuildMaps.c InitMeshInfo.c InitialCondition.c ConvectionRun2d.c
        ConvectionRHS2d.c NcOutput.c)

# set executable
add_executable(Convection2d ${SRC_LIST})

# set static library
target_link_libraries(Convection2d pnetcdf)