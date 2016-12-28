INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(POLYLIB_SOURCE
        polylib.c
        polylib.h)

add_library(polylib STATIC ${POLYLIB_SOURCE})

install(TARGETS polylib ARCHIVE DESTINATION lib)