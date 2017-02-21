INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(PHYSDOMAIN_SOURCE
        pf_phys.c
        pf_phys.h
        pf_fetchBuffer.c
        pf_fetchBuffer.h
        pf_strong_volume_flux2d.c
        pf_strong_volume_flux2d.h
        pf_strong_surface_flux2d.c
        pf_strong_surface_flux2d.h
        pf_add_LDG_solver.c
        pf_add_LDG_solver.h
        pf_strong_viscosity_LDG_flux2d.c
        pf_strong_viscosity_LDG_flux2d.h
        pf_cellMean.c
        pf_cellMean.h
        pf_limiter.c
        pf_limiter.h
        pf_init_file.c
        pf_init_file.h
        pf_openbc.c
        pf_openbc.h)

add_library(physfield STATIC ${PHYSDOMAIN_SOURCE})
add_dependencies(physfield multiregions)
target_link_libraries(physfield multiregions)

install(TARGETS physfield ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)