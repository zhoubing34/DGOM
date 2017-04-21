INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(LIBUTILITIES_SOURCE
        utility.c
        utility.h
        mat_utils.h
        mat_utils.c
        unit_test.h
        unit_test.c
        nc_library.c
        nc_library.h
        arg_section.c
        arg_section.h
        textfile_IO.c)

add_library(utility STATIC ${LIBUTILITIES_SOURCE})
target_link_libraries(utility m lapack blas f2c pnetcdf)# extern lapack library
#if(MKL_FOUND)
#    target_link_libraries(libutility
#        iomp5 pthread m dl mkl_core mkl_intel_lp64 mkl_intel_thread # extern MKL library
#        pnetcdf)
#else(MKL_FOUND)
#    target_link_libraries(libutility
#            m lapack blas f2c pnetcdf)# extern lapack library
#endif(MKL_FOUND)

# install static libraries
install(TARGETS utility ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)