cmake_minimum_required(VERSION 3.1)
project(DGOM)

# set compilers
set(CMAKE_C_FLAGS "-Wall")
set(MPI_PATH $ENV{MPI_DIR})                      # use environmental variable MPI_DIR
set(MPI_C_COMPILER   "${MPI_PATH}/bin/mpicc")    # set mpicc before find mpi package
set(MPI_CXX_COMPILER "${MPI_PATH}/bin/mpicxx")   # set mpicxx before find mpi package

find_package(MPI REQUIRED)
if(NOT MPI_FOUND)
    message(FATAL_ERROR "MPI is not found")
endif()

set(CMAKE_C_COMPILER   ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# build type
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")

# external library
set(THIRD_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/externlib")

# include path
include_directories(
        "${CMAKE_SOURCE_DIR}/library"               # source library
        "${CMAKE_SOURCE_DIR}/include"               # source include
        "${MPI_INCLUDE_PATH}"                       # mpich2 include
        "${THIRD_LIBRARY_PATH}/install/include/"    # external library include
        )

# static library path
LINK_DIRECTORIES("${THIRD_LIBRARY_PATH}/install/lib/")

# set unit test
set(UNITTEST_FOUND 1)
# set precision: double/float
add_definitions(-DDOUBLE_PRECISION)

# add source folder
add_subdirectory(externlib)
#add_subdirectory(src)       # source
add_subdirectory(library)   # static library
add_subdirectory(doc)       # docs
add_subdirectory(test)      # test example

