cmake_minimum_required(VERSION 3.1)
project(DGOM)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

set(CMAKE_C_COMPILER /usr/local/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx)


# 3rdParty library
# Build with original makefile
set(METIS_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/lib/libmetis.a)
set(PARMETIS_LIBRARY ${CMAKE_CURRENT_SOURCE_DIR}/lib/libparmetis.a)

add_custom_target(build_metis ALL
        COMMAND ${CMAKE_MAKE_PROGRAM}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/3rdParty/ParMetis-3.1/METISLib/
        COMMENT "Original metis makefile target")

add_custom_target(build_parmetis ALL
        COMMAND ${CMAKE_MAKE_PROGRAM}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/3rdParty/ParMetis-3.1/ParMETISLib/
        COMMENT "Original parmetis makefile target")

add_library(metis STATIC IMPORTED)
add_library(parmetis STATIC IMPORTED)

set_property(TARGET metis APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(metis PROPERTIES IMPORTED_LOCATION_NOCONFIG "${METIS_LIBRARY}")

set_property(TARGET parmetis APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(parmetis PROPERTIES IMPORTED_LOCATION_NOCONFIG "${PARMETIS_LIBRARY}")

add_dependencies(metis build_metis)
add_dependencies(parmetis build_parmetis)
# source file

add_subdirectory(src)

#set(SRC_LIST
#        src/Utils.c src/LoadBalance2d.c src/FacePair2d.c
#        src/ParallelPairs.c src/BuildMaps2d.c src/StartUp2d.c
#        src/MaxwellsRun2d.c src/MaxwellsMPI2d.c src/MaxwellsDriver2d.c
#        src/MaxwellsRHS2d.c src/InitCPU2d.c)
#
#add_executable(DGOM2 ${SRC_LIST})
#
#include_directories("include")
#target_link_libraries(DGOM2 metis parmetis)