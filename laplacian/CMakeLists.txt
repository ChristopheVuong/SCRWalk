cmake_minimum_required(VERSION 3.0.0)
project(LAPLACIAN)


# add_library(${PROJECT_NAME} src/laplacian.cpp)
add_library(${PROJECT_NAME} src/laplacian_flint.cpp)

#### Flint dependency (for Smith normal form) ####
# SET(CMAKE_FIND_LIBRARY_PREFIXES "lib")
# SET(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
find_library(FLINT flint NAMES flint2 HINTS $ENV{HOME}/flint2)
find_path(FLINT_INCLUDE_DIR flint.h HINTS $ENV{HOME}/flint2)

if (FLINT)
  # include_directories(${FLINT_INCLUDE_DIRS})
  message("Found the library Flint at ${FLINT} !")
  # include_directories(${FLINT})
  target_link_libraries (${PROJECT_NAME} PUBLIC ${FLINT})
  target_include_directories(${PROJECT_NAME} PUBLIC ${FLINT_INCLUDE_DIR})
else()
  message(FATAL_ERROR "No Flint !!!")
endif ()


add_library(sub::liblaplacian ALIAS ${PROJECT_NAME})


# add_library(flint2 STATIC IMPORTED)
# set_property(TARGET flint2 PROPERTY IMPORTED_LOCATION "$ENV{HOME}/flint2")
# include_directories("$ENV{HOME}/flint2")
# target_link_libraries(${PROJECT_NAME}
#     flint
# )

# target_include_directories( ${PROJECT_NAME}
#     PUBLIC ${PROJECT_SOURCE_DIR}/include
# )
