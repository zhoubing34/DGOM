include(ExternalProject)

# CLAPACK-3.1.1.1
set(CLAPACK_PATH  "${THIRD_LIBRARY_PATH}/CLAPACK-3.1.1.1/")
configure_file(${CLAPACK_PATH}/make.inc.example ${CLAPACK_PATH}/make.inc @ONLY)
ExternalProject_Add(
        lapack_lib
        SOURCE_DIR ${CLAPACK_PATH}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND ""
)

# parallel-netcdf-1.7.0
set(PNETCDF_PATH  "${THIRD_LIBRARY_PATH}/parallel-netcdf-1.7.0/")
ExternalProject_Add(
        pnetcdf_lib
        SOURCE_DIR ${PNETCDF_PATH}
        CONFIGURE_COMMAND ./configure --with-mpi=${MPI_PATH} --prefix=${THIRD_LIBRARY_PATH}/install --disable-cxx --disable-fortran MPICC=${CMAKE_C_COMPILER} CC=${CMAKE_C_COMPILER} CFLAGS=-O3
        BUILD_COMMAND make
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND make install
)

# metis-4.0.3
set(METIS_PATH  "${THIRD_LIBRARY_PATH}/parmetis-4.0.3/metis/")
ExternalProject_Add(
        metis_lib
        SOURCE_DIR ${METIS_PATH}
        BINARY_DIR ${METIS_PATH}/build/
        CONFIGURE_COMMAND cmake ../ -DCMAKE_INSTALL_PREFIX=${THIRD_LIBRARY_PATH}/install -DGKLIB_PATH=${METIS_PATH}/GKlib -DSHARED=0 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        BUILD_COMMAND make
        INSTALL_COMMAND make install
)

# parmetis-4.0.3
set(PARMETIS_PATH  "${THIRD_LIBRARY_PATH}/parmetis-4.0.3/")
ExternalProject_Add(
        parmetis_lib
        SOURCE_DIR ${PARMETIS_PATH}
        BINARY_DIR ${PARMETIS_PATH}/build/
        CONFIGURE_COMMAND cmake ../ -DCMAKE_INSTALL_PREFIX=${THIRD_LIBRARY_PATH}/install -DGKLIB_PATH=${METIS_PATH}/GKlib -DMETIS_PATH=${METIS_PATH} -DSHARED=0 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        BUILD_COMMAND make
        INSTALL_COMMAND make install
)