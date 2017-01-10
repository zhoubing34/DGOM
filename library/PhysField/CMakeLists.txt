INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(PHYSDOMAIN_SOURCE
        phys_physField.c
        phys_physField.h
        phys_fetchBuffer.c
        phys_fetchBuffer.h
        phys_strong_volume_flux2d.c
        phys_strong_volume_flux2d.h
        phys_strong_surface_flux2d.c
        phys_strong_surface_flux2d.h
        phys_add_LDG_solver.c
        phys_add_LDG_solver.h
        phys_strong_viscosity_LDG_flux2d.c
        phys_strong_viscosity_LDG_flux2d.h
        phys_cellMean.c
        phys_cellMean.h
        phys_limiter.c
        phys_limiter.h)

add_library(physfield STATIC ${PHYSDOMAIN_SOURCE})
add_dependencies(physfield multiregions)
target_link_libraries(physfield multiregions libutility)

install(TARGETS physfield ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)