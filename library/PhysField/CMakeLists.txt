INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(PHYSDOMAIN_SOURCE
        dg_phys_info.c
        dg_phys_obc.c
        #dg_phys_time_discretize.c
        dg_phys.c
        dg_phys_strong_vol_opt.c
        dg_phys_strong_surf_opt.c
        Limiter/dg_phys_limiter.c
        dg_phys_LDG.c
        dg_phys_strong_LDG_opt.c
        Limiter/dg_phys_limiter_BJ2d.c
        Limiter/dg_phys_indicator_edge.c
        Limiter/dg_phys_indicator_all.c
        )

add_library(phys_field STATIC ${PHYSDOMAIN_SOURCE})
add_dependencies(phys_field multi_area)
target_link_libraries(phys_field multi_area)

install(TARGETS phys_field ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)