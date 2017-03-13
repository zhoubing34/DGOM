# add source file
set(SRC_LIST
        conv_driver.c
        conv_input.c
        conv_init.c
        conv_output.c
        conv_run.c
        conv_rhs.c
        conv_extsol.c
        )

# set executable
add_executable(convection2d ${SRC_LIST})

# set static library
target_link_libraries(convection2d physfield)

# install command
install(TARGETS convection2d RUNTIME DESTINATION bin)