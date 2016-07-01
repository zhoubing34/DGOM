INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/StdRegions)

add_executable(QuadTest QuadTest.c QuadTest.h)
target_link_libraries(QuadTest StdRegions)