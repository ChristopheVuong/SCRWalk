# CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.

cmake_minimum_required(VERSION 3.0.0)
project(Rips-rw)

# Import the cnpy library for saving numpy array
# add_library(cnpy SHARED "cnpy/cnpy.cpp")


# Create the executable
add_executable(${PROJECT_NAME} rips-rw.cpp)
add_executable(${PROJECT_NAME}int rips-rwint.cpp)

#add this library's header folder to the global include set
include_directories("${CMAKE_SOURCE_DIR}/rw/src")
include_directories("${CMAKE_SOURCE_DIR}/laplacian/src")

# Link the static library from subproject1 using its alias sub::lib1
# Link the header only library from subproject2 using its alias sub::lib2
# This will cause the include directories for that target to be added to this project
target_link_libraries(${PROJECT_NAME}
    sub::librw
)
target_link_libraries(${PROJECT_NAME}int
    sub::liblaplacian
)
# target_link_libraries(${PROJECT_NAME} cnpy)
