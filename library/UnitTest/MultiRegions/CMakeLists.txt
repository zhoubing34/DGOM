INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

add_executable(SetCoordTest SetCoordTest.c)

target_link_libraries(SetCoordTest LibUtilities MultiRegions StdRegions)

