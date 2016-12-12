include_directories(${CMAKE_SOURCE_DIR}/library)

set(MULTIREGIONS_TEST_SRC
        MulitRegions_test.c
        MulitRegions_test.h
        VertexSort_test.c
        VertexSort_data.h
        SetTestMultiRegions.c
        SetTestMultiRegions.h
        NodePair_test.h
        NodePair_test.c
        FacePair_test.c
        FacePair_test.h)

add_executable(MulitRegions_Test ${MULTIREGIONS_TEST_SRC})

add_dependencies(MulitRegions_Test MultiRegions)
target_link_libraries(MulitRegions_Test MultiRegions)

install(TARGETS MulitRegions_Test
        RUNTIME DESTINATION test)