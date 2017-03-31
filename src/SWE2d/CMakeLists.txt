set(SWE_LIB
        SWELib/swe_output.c
        SWELib/swe_run.c
        SWELib/swe_rhs.c
        SWELib/swe_flux.c
        SWELib/swe_obc.c
        )
add_library(swe2d_lib STATIC ${SWE_LIB})
add_dependencies(swe2d_lib physfield)
target_link_libraries(swe2d_lib physfield)

set(SRC_LIST
        SWEDriver2d/swe2d.c
        SWEDriver2d/swe_init.c
        SWEDriver2d/swe_input.c
        )

add_executable(swe2d ${SRC_LIST})
target_link_libraries(swe2d swe2d_lib)

# install command
install(TARGETS swe2d RUNTIME DESTINATION bin)
# install examples
install(DIRECTORY example/ DESTINATION test/SWE2d)