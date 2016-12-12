include_directories(${CMAKE_SOURCE_DIR}/library)

set(MULTIREGIONS_TEST_SRC
        MulitRegions_test.c
        MulitRegions_test.h
        VertexSort_test.c
        VertexSort_test.h
        SetTestMultiRegions.c
        SetTestMultiRegions.h
        VarMapPair_test.h
        VarMapPair_test.c
        CellPair_test.c
        CellPair_test.h)

add_executable(MulitRegions_Test ${MULTIREGIONS_TEST_SRC})

add_dependencies(MulitRegions_Test MultiRegions)
target_link_libraries(MulitRegions_Test MultiRegions)

install(TARGETS MulitRegions_Test
        RUNTIME DESTINATION test)