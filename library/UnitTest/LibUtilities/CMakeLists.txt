INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

add_executable(MatrixInverseTest MatrixInverseTest.c)
target_link_libraries(MatrixInverseTest LibUtilities)

#add_executable(MatrixMultiplyTest MatrixMultiplyTest.c)
#target_link_libraries(MatrixMultiplyTest LibUtilities)