INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

# =======================================================================================
# LibUtilities
# =======================================================================================
set(LIBUTILITIES_SOURCE
        LibUtilities/LibUtilities.h
        LibUtilities/Utils.c
        LibUtilities/ParallelPairs.c
        LibUtilities/UTest.c
        LibUtilities/UnstructMesh.c LibUtilities/UnstructMesh.h
        LibUtilities/GenUniformMesh.c LibUtilities/GenUniformMesh.h
        LibUtilities/NetcdfLibrary.c
        LibUtilities/NetcdfLibrary.h
        LibUtilities/SlopeLimiter.c LibUtilities/SlopeLimiter.h)
add_library(LibUtilities STATIC ${LIBUTILITIES_SOURCE})
TARGET_LINK_LIBRARIES(LibUtilities
        iomp5 pthread m dl mkl_core mkl_intel_lp64 mkl_intel_thread # extern MKL library
        pnetcdf)

# =======================================================================================
# polylib
# =======================================================================================
set(POLYLIB_SOURCE
        Polylib/polylib.c
        Polylib/polylib.h)

add_library(Polylib STATIC ${POLYLIB_SOURCE})


# =======================================================================================
# StdRegions
# =======================================================================================

set(STDREGIONS_SOURCE
        StdRegions/StdRegions.h
        StdRegions/Triangle.c
        StdRegions/Quadrilateral.c)

add_library(StdRegions STATIC ${STDREGIONS_SOURCE})
add_dependencies(StdRegions LibUtilities Polylib)
target_link_libraries(StdRegions LibUtilities Polylib)

# =======================================================================================
# LocalRegions
# =======================================================================================
set(LOCAL_REGIONS_SOURCE
        LocalRegions/LocalRegions.h
        LocalRegions/Mapping.c)

add_library(LocalRegions STATIC ${LOCAL_REGIONS_SOURCE})
add_dependencies(LocalRegions LibUtilities Polylib)
target_link_libraries(LocalRegions LibUtilities Polylib)

# =======================================================================================
# MultiRegions
# =======================================================================================
set(MULTI_REGIONS_SOURCE
        MultiRegions/MultiRegions.c MultiRegions/MultiRegions.h
        MultiRegions/SetFacePair.c
        MultiRegions/SetNodePair.c
        MultiRegions/SetElemetPair.c
        MultiRegions/LoadBalance.c)

add_library(MultiRegions STATIC ${MULTI_REGIONS_SOURCE})
add_dependencies(MultiRegions LocalRegions)
target_link_libraries(MultiRegions LibUtilities
        LocalRegions StdRegions
        parmetis metis)

# =======================================================================================
# MultiRegions
# =======================================================================================
set(PHYSDOMAIN_SOURCE
        PhysDomain/PhysDomain.c
        PhysDomain/PhysDomain.h)

add_library(PhysDomain STATIC ${PHYSDOMAIN_SOURCE})
add_dependencies(PhysDomain MultiRegions)
target_link_libraries(PhysDomain MultiRegions LibUtilities)

# =======================================================================================
# UnitTest, configure with `-DUnitTest=1`
# =======================================================================================
set(UnitTest 1) # turn on the model

if(DEFINED UnitTest)
    add_subdirectory(UnitTest)
endif()

# =======================================================================================
# install command
# =======================================================================================
install(TARGETS
        LibUtilities Polylib StdRegions
        LocalRegions MultiRegions PhysDomain
        ARCHIVE DESTINATION lib)

