# add source file
set(SRC_LIST
        conv_driver2d.c
        conv_driver2d.h
        conv_mesh.c
        conv_mesh.h
        conv_getparameter.c
        conv_getparameter.h
        conv_intilization.c
        conv_intilization.h
        conv_output.c
        conv_output.h
        conv_run.c
        conv_run.h
        conv_rhs.c
        conv_rhs.h
        conv_extsol.c
        conv_extsol.h)

# set executable
add_executable(convection2d ${SRC_LIST})

# set static library
target_link_libraries(convection2d pnetcdf physfield)

# install command
install(TARGETS convection2d RUNTIME DESTINATION bin)