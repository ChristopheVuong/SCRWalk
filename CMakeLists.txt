cmake_minimum_required(VERSION 3.0.0)
project(SCRWalk)

#### GUDHI requires Boost - Alpha complex requires Boost program_options thread ####
find_package(Boost 1.48.0 REQUIRED COMPONENTS program_options thread)

# BOOST ISSUE result_of vs C++11
add_definitions(-DBOOST_RESULT_OF_USE_DECLTYPE)
# BOOST ISSUE with Libraries name resolution under Windows
add_definitions(-DBOOST_ALL_NO_LIB)
# problem with Visual Studio link on Boost program_options
add_definitions( -DBOOST_ALL_DYN_LINK )
# problem on Mac with boost_system and boost_thread
add_definitions( -DBOOST_SYSTEM_NO_DEPRECATED )

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

#### GUDHI Alpha complex requires CGAL >= 4.7.0 ####
find_package(CGAL QUIET)

# Requires CGAL versions > 4.11
if (NOT CGAL_FOUND)
  message(FATAL_ERROR "++ CGAL is required for this example.")
endif()

message("++ CGAL Version is ${CGAL_VERSION}")

if (CGAL_VERSION VERSION_LESS 4.11.0)
  message(FATAL_ERROR "++ Your CGAL is considered too old to be used by Gudhi. Must be > 4.11.")
endif()

include( ${CGAL_USE_FILE} )

#### Optional GMP and GMPXX for CGAL ####
find_package(GMP)
if(GMP_FOUND)
  include_directories(${GMP_INCLUDE_DIR})
  find_package(GMPXX)
  if(GMPXX_FOUND)
    include_directories(${GMPXX_INCLUDE_DIR})
  endif()
endif()

#### GUDHI Alpha complex requires Eigen3 ####
find_package(Eigen3 3.1.0 REQUIRED)
if (EIGEN3_FOUND)
  include( ${EIGEN3_USE_FILE} )
endif (EIGEN3_FOUND)

#### Optional TBB for CGAL and GUDHI ####
set(TBB_FIND_QUIETLY ON)
find_package(TBB)
if (TBB_FOUND)
  include(${TBB_USE_FILE})
endif()

message("GUDHI_INCLUDE_DIRS = ${GUDHI_INCLUDE_DIRS}")

#### GUDHI Alpha complex requires GUDHI >= 2.0.0 ####
find_package(GUDHI 2.0.0 REQUIRED)
message("GUDHI_VERSION = ${GUDHI_VERSION}")
message("GUDHI_INCLUDE_DIRS = ${GUDHI_INCLUDE_DIRS}")
include_directories(${GUDHI_INCLUDE_DIRS})


# Recurse into the subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
# add_subdirectory (src)
# add_subdirectory (tests)

# Add a library with the above sources
# Add sub directories
add_subdirectory(rw)
add_subdirectory(laplacian)
# add subdirectory binary executable
add_subdirectory(examples)

# add_subdirectory(cnpy)
# add_library(${PROJECT_NAME} rw/src/rwZ2chains.cpp rw/src/rwZchains.cpp)
# add_library(sub::librw ALIAS ${PROJECT_NAME})



include(CTest)
enable_testing()

