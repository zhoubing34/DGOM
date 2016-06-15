INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)
#set(THIRD_LIBRARY_PATH "${CMAKE_SOURCE_DIR}/3rdParty/install/lib/")
#LINK_DIRECTORIES(${THIRD_LIBRARY_PATH})

# LibUtilities
set(LIBUTILITIES_HEADER
        LibUtilities/LibUtilities.h)

set(LIBUTILITIES_SOURCE LibUtilities/Utils.c
        LibUtilities/ParallelPairs.c)
add_library(LibUtilities STATIC ${LIBUTILITIES_SOURCE})
TARGET_LINK_LIBRARIES(LibUtilities
        iomp5 pthread m dl mkl_core mkl_intel_lp64 mkl_intel_thread) # extern MKL library

# polylib
set(POLYLIB_HEADER Polylib/polylib.h)
set(POLYLIB_SOURCE Polylib/polylib.c Polylib/polylib.h)

add_library(Polylib STATIC ${POLYLIB_SOURCE})

# StdRegions
set(STDREGIONS_HEADER
        StdRegions/StdRegions.h)

set(STDREGIONS_SOURCE
        StdRegions/Triangle.c
        StdRegions/Quadrilateral.c)

add_library(StdRegions STATIC ${STDREGIONS_SOURCE})
#add_dependencies(StdRegions LibUtilities Polylib)
##TARGET_LINK_LIBRARIES(StdRegions blas f2c lapack) # extern library -lblas -lf2c -lapack
target_link_libraries(StdRegions LibUtilities Polylib)
#
## LocalRegions
#set(LOCAL_REGIONS_SOURCE LocalRegions/TriMapping.c)
#
#add_library(LocalRegions STATIC ${LOCAL_REGIONS_SOURCE})
#add_dependencies(LocalRegions LibUtilities Polylib)
#target_link_libraries(LocalRegions Polylib)


#
