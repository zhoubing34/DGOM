set(SRC_LIST
        ConvectionDriver2d.c Mesh2d.c ../Utilities/Utils.c ../Utilities/ParallelPairs.c
        PairFace.c SetUp.c BuildMaps.c InitMeshInfo.c InitialCondition.c ConvectionRun2d.c
        ConvectionRHS2d.c)

add_executable(Convection2d ${SRC_LIST})

include_directories("../../include" "/usr/include")
target_link_libraries(Convection2d metis parmetis)