include_directories(${CMAKE_SOURCE_DIR}/library)

set(PHYSDOMAIN_TEST_SRC
        PhysDomain_test.c
        PhysDomain_test.h)

add_executable(PhysDomain_Test ${PHYSDOMAIN_TEST_SRC})
add_dependencies(PhysDomain_Test PhysDomain)
target_link_libraries(PhysDomain_Test PhysDomain)

install(TARGETS PhysDomain_Test
        RUNTIME DESTINATION test)
