set(SRC_LIST SWE2d_Test.c)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

add_executable(SWE2d_Test ${SRC_LIST})

target_link_libraries(SWE2d_Test pnetcdf)
target_link_libraries(SWE2d_Test parmetis)
target_link_libraries(SWE2d_Test metis)
target_link_libraries(SWE2d_Test PhysDomain)