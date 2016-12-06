INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(LIBUTILITIES_SOURCE
        LibUtilities.h
        Utils.c
        ParallelPairs.c
        UTest.c
        UnstructMesh.c
        UnstructMesh.h
        GenUniformMesh.c
        GenUniformMesh.h
        NetcdfLibrary.c
        NetcdfLibrary.h
        SlopeLimiter.c
        SlopeLimiter.h)

add_library(LibUtilities STATIC ${LIBUTILITIES_SOURCE})

TARGET_LINK_LIBRARIES(LibUtilities
        iomp5 pthread m dl mkl_core mkl_intel_lp64 mkl_intel_thread # extern MKL library
        pnetcdf)

install(TARGETS LibUtilities ARCHIVE DESTINATION lib)