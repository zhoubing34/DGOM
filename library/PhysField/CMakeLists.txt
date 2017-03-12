INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(PHYSDOMAIN_SOURCE
        dg_phys.c
        dg_phys_strong_vol_opt.c
        dg_phys_strong_surf_opt.c
        #pf_add_LDG_solver.c
        #pf_strong_viscosity_LDG_flux2d.c
        #pf_cellMean.c
        #pf_init_file.c
        #pf_openbc.c
        #Limiter/pf_limit_BJ2d.c)
        )

add_library(physfield STATIC ${PHYSDOMAIN_SOURCE})
add_dependencies(physfield multiregions)
target_link_libraries(physfield multiregions)

install(TARGETS physfield ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)