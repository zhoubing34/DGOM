add_subdirectory(Polylib)
add_subdirectory(LibUtilities)
add_subdirectory(StdRegions)
add_subdirectory(LocalRegions)
add_subdirectory(MultiRegions)
add_subdirectory(PhysDomain)

# =======================================================================================
# UnitTest, configure with `-DUnitTest=1`
# =======================================================================================
set(UnitTest 1) # turn on the model

if(DEFINED UnitTest)
    add_subdirectory(UnitTest)
endif()