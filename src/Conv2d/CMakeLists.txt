set(CONV_LIB
        ConvLib/conv_output.c
        ConvLib/conv_run.c
        ConvLib/conv_rhs.c
        )

add_library(conv2d_lib STATIC ${CONV_LIB})
add_dependencies(conv2d_lib phys_field)
target_link_libraries(conv2d_lib phys_field)

# standard conv2d solver
set(SRC_LIST
        Conv2d_st/conv2d_st.c
        Conv2d_st/conv_input_st.c
        Conv2d_st/conv_init_st.c
        Conv2d_st/conv_extsol.c
        )
add_executable(conv2d_st ${SRC_LIST}) # set executable
target_link_libraries(conv2d_st conv2d_lib) # set static library
install(TARGETS conv2d_st RUNTIME DESTINATION bin) # install command

# user specific conv2d solver
set(SRC_LIST
        ConvDriver2d/conv2d.c
        ConvDriver2d/conv_input.c
        ConvDriver2d/conv_init.c
        )

add_executable(conv2d ${SRC_LIST}) # set executable
target_link_libraries(conv2d conv2d_lib) # set static library
install(TARGETS conv2d RUNTIME DESTINATION bin) # install command

install(DIRECTORY example/ DESTINATION test/Conv2d)