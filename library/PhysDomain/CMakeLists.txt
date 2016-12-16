INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(PHYSDOMAIN_SOURCE
        PhysDomain.c
        PhysDomain.h)
#        SlopeLimiter.c
#        SlopeLimiter.h)

add_library(PhysDomain STATIC ${PHYSDOMAIN_SOURCE})
add_dependencies(PhysDomain MultiRegions)
target_link_libraries(PhysDomain MultiRegions LibUtilities)

install(TARGETS PhysDomain ARCHIVE DESTINATION lib)

if(UNITTEST_FOUND)
    add_subdirectory(test)
endif(UNITTEST_FOUND)