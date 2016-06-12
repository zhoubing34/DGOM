# add source file
set(SRC_LIST
        Mesh2d.c ../../library/LibUtilities/Utils.c LoadBalance2d.c FacePair2d.c
        ../../library/LibUtilities/ParallelPairs.c BuildMaps2d.c StartUp2d.c
        MaxwellsRun2d.c MaxwellsMPI2d.c MaxwellsDriver2d.c
        MaxwellsRHS2d.c InitCPU2d.c)

# set executable
add_executable(Maxwells2d ${SRC_LIST})

# add static library
target_link_libraries(Maxwells2d parmetis)
target_link_libraries(Maxwells2d metis)