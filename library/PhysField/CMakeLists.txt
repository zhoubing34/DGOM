INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(PHYSDOMAIN_SOURCE
        phs_physField.c
        phs_physField.h
        phys_fetchBuffer.c
        phys_fetchBuffer.h)

add_library(PhysField STATIC ${PHYSDOMAIN_SOURCE})
add_dependencies(PhysField MultiRegions)
target_link_libraries(PhysField MultiRegions LibUtilities)

install(TARGETS PhysField ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)