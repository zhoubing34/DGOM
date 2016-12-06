INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/library)

set(POLYLIB_SOURCE
        polylib.c
        polylib.h)

add_library(Polylib STATIC ${POLYLIB_SOURCE})

install(TARGETS Polylib ARCHIVE DESTINATION lib)