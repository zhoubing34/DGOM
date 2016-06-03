cmake_minimum_required(VERSION 3.1)
project(DGOM)

#set constant
set(PNETCDF_PATH "${CMAKE_SOURCE_DIR}/3rdParty/install")  # pnetcdf
set(MPI_PATH "/usr/local/Cellar/mpich2/1.4.1/")     	# mpich2
set(METIS_PATH "${CMAKE_SOURCE_DIR}/3rdParty/install")	# metis and parmetis

# set compilers
set(CMAKE_C_FLAGS "-Wall")
set(CMAKE_C_COMPILER /usr/local/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx)

# add header path
include_directories("${CMAKE_SOURCE_DIR}/include"   # source include
        "${MPI_PATH}/include"                       # mpi include
        "${PNETCDF_PATH}/include/")                 # third library include

# add static library path
LINK_DIRECTORIES("${NETCDF_PATH}/lib/" "${METIS_PATH}/lib")

# source file
add_subdirectory(src)

# add docs
add_subdirectory(doc)
