INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

#
add_executable(SetCoordTest SetCoordTest.c)
target_link_libraries(SetCoordTest LibUtilities MultiRegions StdRegions)

#
add_executable(LoadBalanceTest LoadBalanceTest.c)
target_link_libraries(LoadBalanceTest LibUtilities MultiRegions)


add_executable(SetFacePairTest SetFacePairTest.c)
target_link_libraries(SetFacePairTest LibUtilities MultiRegions)


add_executable(SetNodePairTest SetNodePairTest.c)
target_link_libraries(SetNodePairTest LibUtilities MultiRegions)

add_executable(SetVolumeGeoTest SetVolumeGeoTest.c)
target_link_libraries(SetVolumeGeoTest LibUtilities MultiRegions)

