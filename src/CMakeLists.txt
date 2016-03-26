set(SRC_LIST
        Mesh2d.c Utils.c LoadBalance2d.c FacePair2d.c
        ParallelPairs.c BuildMaps2d.c StartUp2d.c
        MaxwellsRun2d.c MaxwellsMPI2d.c MaxwellsDriver2d.c
        MaxwellsRHS2d.c InitCPU2d.c)

add_executable(DGOM2 ${SRC_LIST})

include_directories("../include")
target_link_libraries(DGOM2 metis parmetis)