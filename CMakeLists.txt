cmake_minimum_required(VERSION 3.0.0)
project(SCRWalk)


# Recurse into the subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
# add_subdirectory (src)
# add_subdirectory (tests)

# Add a library with the above sources
# Add sub directories
add_subdirectory(rw)
# add_subdirectory(laplacian)
# add subdirectory binary executable
add_subdirectory(examples)
# add_library(${PROJECT_NAME} rw/src/rwZ2chains.cpp rw/src/rwZchains.cpp)
# add_library(sub::librw ALIAS ${PROJECT_NAME})



include(CTest)
enable_testing()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
