INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/PhysDomain)

#
add_executable(SetMapOUT2dTest SetMapOUT2dTest.c)
target_link_libraries(SetMapOUT2dTest PhysDomain)