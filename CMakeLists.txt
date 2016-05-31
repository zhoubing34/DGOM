cmake_minimum_required(VERSION 3.1)
project(DGOM)

#set constant
set(NETCDF_PATH "/usr/local/Cellar/netcdf/3.6.2/")  # netcdf
set(MPI_PATH "/usr/local/Cellar/mpich2/1.4.1/")     # mpich2
set(METIS_PATH "${CMAKE_SOURCE_DIR}")               # metis and parmetis

# set compilers
set(CMAKE_C_FLAGS "-Wall")
set(CMAKE_C_COMPILER /usr/local/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx)

# add header path
include_directories("${CMAKE_SOURCE_DIR}/include"
        "${MPI_PATH}/include"
        "${NETCDF_PATH}/include/")

# add static library path
LINK_DIRECTORIES("${NETCDF_PATH}/lib/" "${METIS_PATH}/lib")

# source file
add_subdirectory(src)
