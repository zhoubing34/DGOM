INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

# LibUtilities
set(LIBUTILITIES_HEADER
        LibUtilities/LibUtilities.h)

set(LIBUTILITIES_SOURCE
        LibUtilities/ParallelPairs.c
        LibUtilities/Utils.c)

add_library(LibUtilities STATIC ${LIBUTILITIES_SOURCE})

# polylib
set(POLYLIB_HEADER
        Polylib/polylib.h)

set(POLYLIB_SOURCE
        Polylib/polylib.c
        Polylib/polylib.h)
add_library(Polylib STATIC ${POLYLIB_SOURCE})

# StdRegions
set(STDREGIONS_HEADER
        StdRegions/StdRegions.h)

set(STDREGIONS_SOURCE
        StdRegions/Triangle.c
        StdRegions/Quadrilateral.c)

add_library(StdRegions STATIC ${STDREGIONS_SOURCE})

TARGET_LINK_LIBRARIES(StdRegions LibUtilities Polylib)

# LocalRegions
set(LOCAL_REGIONS_SOURCE LocalRegions/TriMapping.c)

add_library(LocalRegions STATIC ${LOCAL_REGIONS_SOURCE})
TARGET_LINK_LIBRARIES(LocalRegions Polylib)

#
