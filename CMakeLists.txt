project(mpm LANGUAGES CXX)

# Require C++14-compliant compiler; only available for CMake v. 3.1 and up
set(CMAKE_CXX_STANDARD 14)

cmake_minimum_required(VERSION 3.1)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

# General compile settings
IF (NOT CMAKE_BUILD_TYPE)
 SET(CMAKE_BUILD_TYPE "Debug")
 #SET(CMAKE_BUILD_TYPE "Release")
ENDIF (NOT CMAKE_BUILD_TYPE)

# GNU Specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
endif()

# CMake Modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Boost Archive
find_package(Boost REQUIRED COMPONENTS system)
include_directories(${BOOST_INCLUDE_DIRS})
link_libraries(${Boost_SYSTEM_LIBRARY})

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# TBB
find_package(TBB REQUIRED)
add_definitions(${TBB_DEFINITIONS})
include_directories(${TBB_INCLUDE_DIRS})
link_libraries(${TBB_LIBRARIES})

# pthreads
find_package (Threads)
link_libraries (${CMAKE_THREAD_LIBS_INIT})

# Include directories
include_directories(
  ${mpm_SOURCE_DIR}/include/
  ${mpm_SOURCE_DIR}/external/
  ${mpm_SOURCE_DIR}/tests/include/
)

# mpm executable
SET(mpm_src
  ${mpm_SOURCE_DIR}/src/affine_transform.cc
  ${mpm_SOURCE_DIR}/src/cell.cc
  ${mpm_SOURCE_DIR}/src/io.cc
  ${mpm_SOURCE_DIR}/src/logger.cc
  ${mpm_SOURCE_DIR}/src/material.cc
  ${mpm_SOURCE_DIR}/src/mpm.cc
  ${mpm_SOURCE_DIR}/src/node.cc
  ${mpm_SOURCE_DIR}/src/particle.cc
  ${mpm_SOURCE_DIR}/src/read_mesh.cc
  ${mpm_SOURCE_DIR}/src/shapefn.cc
  ${mpm_SOURCE_DIR}/src/main.cc
)
add_executable(mpm
  ${mpm_src}
)

# Unit test
SET(test_src
  ${mpm_SOURCE_DIR}/src/affine_transform.cc
  ${mpm_SOURCE_DIR}/src/cell.cc
  ${mpm_SOURCE_DIR}/src/io.cc
  ${mpm_SOURCE_DIR}/src/logger.cc
  ${mpm_SOURCE_DIR}/src/material.cc
  ${mpm_SOURCE_DIR}/src/mpm.cc
  ${mpm_SOURCE_DIR}/src/node.cc
  ${mpm_SOURCE_DIR}/src/particle.cc
  ${mpm_SOURCE_DIR}/src/read_mesh.cc
  ${mpm_SOURCE_DIR}/src/shapefn.cc
  ${mpm_SOURCE_DIR}/tests/cell_container_test.cc
  ${mpm_SOURCE_DIR}/tests/cell_test.cc
  ${mpm_SOURCE_DIR}/tests/hex_quadrature_test.cc  
  ${mpm_SOURCE_DIR}/tests/hex_shapefn_test.cc
  ${mpm_SOURCE_DIR}/tests/io_test.cc
  ${mpm_SOURCE_DIR}/tests/material/linear_elastic_test.cc
  ${mpm_SOURCE_DIR}/tests/mesh_test.cc
  #${mpm_SOURCE_DIR}/tests/mpm_explicit_test.cc
  #${mpm_SOURCE_DIR}/tests/mpm_explicit_test_unitcell.cc
  ${mpm_SOURCE_DIR}/tests/node_container_test.cc
  ${mpm_SOURCE_DIR}/tests/node_handler_test.cc
  ${mpm_SOURCE_DIR}/tests/node_test.cc
  ${mpm_SOURCE_DIR}/tests/particle_container_test.cc
  ${mpm_SOURCE_DIR}/tests/particle_test.cc
  ${mpm_SOURCE_DIR}/tests/quad_quadrature_test.cc    
  ${mpm_SOURCE_DIR}/tests/quad_shapefn_test.cc
  ${mpm_SOURCE_DIR}/tests/read_mesh_ascii_test.cc
  ${mpm_SOURCE_DIR}/tests/test_main.cc
)   
add_executable(mpmtest
  ${test_src}
)
add_test(NAME mpmtest COMMAND $<TARGET_FILE:mpmtest>)
enable_testing()

# Coverage
find_package(codecov)
add_coverage(mpm)
add_coverage(mpmtest)
