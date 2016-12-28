INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(LIBUTILITIES_SOURCE
        LibUtilities.h
        Utils.h
        Utils.c
        UTest.h
        UTest.c
        nc_library.c
        nc_library.h)

add_library(libutility STATIC ${LIBUTILITIES_SOURCE})

TARGET_LINK_LIBRARIES(libutility
        iomp5 pthread m dl mkl_core mkl_intel_lp64 mkl_intel_thread # extern MKL library
        pnetcdf)

# install static libraries
install(TARGETS libutility ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)