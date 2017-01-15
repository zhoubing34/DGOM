set(SRC_LIST
        swe_driver2d.c
        swe_driver2d.h
        swe_init.c
        swe_init.h
        swe_input.c
        swe_input.h
        swe_read_command.c
        swe_read_command.h
        swe_mesh.c
        swe_mesh.h
        swe_output.c
        swe_output.h
        swe_run.c
        swe_run.h
        swe_rhs.c
        swe_rhs.h
        )

add_executable(swe2d ${SRC_LIST})

target_link_libraries(swe2d pnetcdf)
target_link_libraries(swe2d parmetis)
target_link_libraries(swe2d metis)
target_link_libraries(swe2d physfield)

# install command
install(TARGETS swe2d RUNTIME DESTINATION bin)