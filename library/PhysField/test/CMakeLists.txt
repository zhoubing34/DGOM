include_directories(${CMAKE_SOURCE_DIR}/library)

set(PHYSDOMAIN_TEST_SRC
        phys_test.c
        phys_test.h
        phys_nodeFetch_test.c
        phys_nodeFetch_test.h
        phys_strong_volume_flux2d_test.c
        phys_strong_volume_flux2d_test.h
        phys_strong_surface_flux2d_test.c
        phys_strong_surface_flux2d_test.h)

add_executable(PhysField_Test ${PHYSDOMAIN_TEST_SRC})
add_dependencies(PhysField_Test PhysField)
target_link_libraries(PhysField_Test PhysField)

install(TARGETS PhysField_Test
        RUNTIME DESTINATION test)
