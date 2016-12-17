include_directories(${CMAKE_SOURCE_DIR}/library)

set(STDREGIONS_SOURCE
        sc_stdcell.h
        sc_stdcell.c
        sc_triangle.c
        sc_quadrilateral.c)

add_library(StandCell STATIC ${STDREGIONS_SOURCE})
add_dependencies(StandCell LibUtilities Polylib)
target_link_libraries(StandCell LibUtilities Polylib)

install(TARGETS StandCell ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)