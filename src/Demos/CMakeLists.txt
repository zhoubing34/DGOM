# set source file
set(DEMOS_SOURCE test.c)

include_directories(${CMAKE_SOURCE_DIR}/library)

# add executable file
add_executable(Demos ${DEMOS_SOURCE})

target_link_libraries(Demos LibUtilities StdRegions)
#        "/opt/intel/composer_xe_2015.5.222/mkl/lib/libmkl_core.a"
#        "/opt/intel/composer_xe_2015.5.222/mkl/lib/libmkl_intel_lp64.a"
#        "/opt/intel/composer_xe_2015.5.222/mkl/lib/libmkl_intel_thread.a")