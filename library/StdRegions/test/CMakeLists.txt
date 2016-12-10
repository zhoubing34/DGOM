INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library
        ${CMAKE_SOURCE_DIR}/StdRegions)

set(STDREGIONS_SCR
        StdRegions_Test.c
        StdRegions_Test.h)

add_executable(StdRegions_Test ${STDREGIONS_SCR})

add_dependencies(StdRegions_Test LibUtilities StdRegions)
target_link_libraries(StdRegions_Test LibUtilities StdRegions)

install(TARGETS StdRegions_Test
        RUNTIME DESTINATION test)
