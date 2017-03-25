include_directories(${CMAKE_SOURCE_DIR}/library)

set(PHYSDOMAIN_TEST_SRC
        dg_phys_test.c
        dg_phys_strong_vol_opt_test.c
        dg_phys_strong_surf_opt_test.c
        #pf_strong_viscosity_LDG_flux2d_test.c
        #pf_cellMean_test.c
        dg_phys_limiter_test.c
        dg_phys_obc_test.c
        )

add_executable(physfield_test ${PHYSDOMAIN_TEST_SRC})
add_dependencies(physfield_test physfield)
target_link_libraries(physfield_test physfield)

install(TARGETS physfield_test
        RUNTIME DESTINATION test)
