cmake_minimum_required(VERSION 3.1)
project(DGOM)

#set constant
set(PNETCDF_PATH "${CMAKE_SOURCE_DIR}/3rdParty/install")  # pnetcdf
set(METIS_PATH "${CMAKE_SOURCE_DIR}/3rdParty/install")	# metis and parmetis

# set mpi compilers
set(CMAKE_C_FLAGS "-Wall")
set(MPI_PATH $ENV{MPI_DIR})                      # use environmental variable MPI_DIR
set(MPI_C_COMPILER "${MPI_PATH}/bin/mpicc")      # set mpicc before find mpi package
set(MPI_CXX_COMPILER "${MPI_PATH}/bin/mpicxx")   # set mpicxx before find mpi package

find_package(MPI REQUIRED)
if(NOT MPI_FOUND)
    message(FATAL_ERROR "MPI is not found")
endif()

set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# add header path
include_directories("${CMAKE_SOURCE_DIR}/include"   # source include
        "${MPI_INCLUDE_PATH}"                       # mpi include
        "${PNETCDF_PATH}/include/")                 # third library include

# add static library path
LINK_DIRECTORIES("${NETCDF_PATH}/lib/" "${METIS_PATH}/lib")

# source file
SUBDIRS(src library)

# add docs
add_subdirectory(doc)
