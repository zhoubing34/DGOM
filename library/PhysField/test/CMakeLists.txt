include_directories(${CMAKE_SOURCE_DIR}/library)

set(PHYSDOMAIN_TEST_SRC
        dg_phys_test.c
        dg_phys_strong_vol_opt_test.c
        dg_phys_strong_surf_opt_test.c
        dg_phys_strong_LDG_opt_test.c
        dg_phys_limiter_test.c
        dg_phys_obc_test.c
        )

add_executable(phys_field_test ${PHYSDOMAIN_TEST_SRC})
add_dependencies(phys_field_test phys_field)
target_link_libraries(phys_field_test phys_field)

install(TARGETS phys_field_test RUNTIME DESTINATION test)
