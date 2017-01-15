include_directories(${CMAKE_SOURCE_DIR}/library)

set(PHYSDOMAIN_TEST_SRC
        pf_test.c
        pf_test.h
        pf_nodeFetch_test.c
        pf_nodeFetch_test.h
        pf_strong_volume_flux2d_test.c
        pf_strong_volume_flux2d_test.h
        pf_strong_surface_flux2d_test.c
        pf_strong_surface_flux2d_test.h
        pf_strong_viscosity_LDG_flux2d_test.h
        pf_strong_viscosity_LDG_flux2d_test.c
        pf_cellMean_test.c
        pf_cellMean_test.h
        pf_cellFetch_test.c
        pf_cellFetch_test.h
        pf_limiter_test.c
        pf_limiter_test.h)

add_executable(physfield_test ${PHYSDOMAIN_TEST_SRC})
add_dependencies(physfield_test physfield)
target_link_libraries(physfield_test physfield)

install(TARGETS physfield_test
        RUNTIME DESTINATION test)
