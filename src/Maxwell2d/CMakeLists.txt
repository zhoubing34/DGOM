set(SRC_LIST
        Mesh2d.c ../Utilities/Utils.c LoadBalance2d.c FacePair2d.c
        ../Utilities/ParallelPairs.c BuildMaps2d.c StartUp2d.c
        MaxwellsRun2d.c MaxwellsMPI2d.c MaxwellsDriver2d.c
        MaxwellsRHS2d.c InitCPU2d.c)

add_executable(Maxwells2d ${SRC_LIST})

include_directories("../../include")
target_link_libraries(Maxwells2d metis parmetis)