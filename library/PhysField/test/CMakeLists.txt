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

add_executable(physfield_test ${PHYSDOMAIN_TEST_SRC})
add_dependencies(physfield_test physfield)
target_link_libraries(physfield_test physfield)

install(TARGETS physfield_test
        RUNTIME DESTINATION test)
