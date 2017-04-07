set(CONV_LIB
        ConvLib/conv_output.c
        ConvLib/conv_run.c
        ConvLib/conv_rhs.c
        )

add_library(conv2d_lib STATIC ${CONV_LIB})
add_dependencies(conv2d_lib physfield)
target_link_libraries(conv2d_lib physfield)

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