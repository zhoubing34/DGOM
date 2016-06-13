# set source file
set(DEMOS_SOURCE test.c)

include_directories(${CMAKE_SOURCE_DIR}/library)

# add executable file
add_executable(Demos ${DEMOS_SOURCE})

target_link_libraries(Demos LibUtilities StdRegions)