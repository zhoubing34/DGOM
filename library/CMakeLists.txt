INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

# =======================================================================================
# LibUtilities
# =======================================================================================
set(LIBUTILITIES_HEADER
        LibUtilities/LibUtilities.h)

set(LIBUTILITIES_SOURCE
        LibUtilities/Utils.c
        LibUtilities/ParallelPairs.c
        LibUtilities/UTest.c)
add_library(LibUtilities STATIC ${LIBUTILITIES_SOURCE})
TARGET_LINK_LIBRARIES(LibUtilities
        iomp5 pthread m dl mkl_core mkl_intel_lp64 mkl_intel_thread) # extern MKL library

# =======================================================================================
# polylib
# =======================================================================================
set(POLYLIB_HEADER Polylib/polylib.h)
set(POLYLIB_SOURCE Polylib/polylib.c Polylib/polylib.h)

add_library(Polylib STATIC ${POLYLIB_SOURCE})


# =======================================================================================
# StdRegions
# =======================================================================================
set(STDREGIONS_HEADER
        StdRegions/StdRegions.h)

set(STDREGIONS_SOURCE
        StdRegions/Triangle.c
        StdRegions/Quadrilateral.c)

add_library(StdRegions STATIC ${STDREGIONS_SOURCE})
add_dependencies(StdRegions LibUtilities Polylib)
target_link_libraries(StdRegions LibUtilities Polylib)

# =======================================================================================
# LocalRegions
# =======================================================================================
set(LOCAL_REGIONS_SOURCE LocalRegions/Mapping.c)

add_library(LocalRegions STATIC ${LOCAL_REGIONS_SOURCE})
add_dependencies(LocalRegions LibUtilities Polylib)
target_link_libraries(LocalRegions LibUtilities Polylib)

# =======================================================================================
# MultiRegions
# =======================================================================================
set(MULTI_REGIONS_SOURCE
        MultiRegions/MultiRegions.c
        MultiRegions/MultiRegions.h
        MultiRegions/SetFacePair.c
        MultiRegions/SetNodePair.c
        MultiRegions/LoadBalance.c)

add_library(MultiRegions STATIC ${MULTI_REGIONS_SOURCE})
add_dependencies(MultiRegions LocalRegions)
target_link_libraries(MultiRegions LibUtilities
        LocalRegions StdRegions
        parmetis metis)

# =======================================================================================
# UnitTest, configure with `-DUnitTest=1`
# =======================================================================================
#set(UnitTest 1) # turn on the model

if(DEFINED UnitTest)
    add_subdirectory(UnitTest)
endif()

