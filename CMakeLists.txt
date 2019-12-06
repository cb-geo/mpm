project(mpm LANGUAGES CXX)

# Require C++14-compliant compiler; only available for CMake v3.10 and up
set(CMAKE_CXX_STANDARD 14)

cmake_minimum_required(VERSION 3.10)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

# General compile settings
IF (NOT CMAKE_BUILD_TYPE)
 #SET(CMAKE_BUILD_TYPE "Debug")
 SET(CMAKE_BUILD_TYPE "Release")
ENDIF (NOT CMAKE_BUILD_TYPE)

# GNU specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
endif()

# Intel specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()

# Clang specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template")
endif()

# CMake seems to have no way to enable/disable testing per subproject,
# so we provide an option similar to BUILD_TESTING, but just for MPM.
option(MPM_BUILD_TESTING "enable testing for mpm" ON)

# CMake option for building mpm as a shared library for CI
option(MPM_BUILD_LIB "enable libmpm" OFF)

# CMake Modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Boost Archive
find_package(Boost REQUIRED COMPONENTS filesystem system)
include_directories(${BOOST_INCLUDE_DIRS})
link_libraries(${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# MPI
find_package(MPI)
if (MPI_FOUND)
  add_definitions("-DUSE_MPI")
  include_directories(${MPI_CXX_INCLUDE_DIRS})
  link_libraries(${MPI_CXX_LIBRARIES})
endif()

# HDF5
ENABLE_LANGUAGE(C)
find_package(HDF5 COMPONENTS CXX HL)
if (HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
  link_libraries(${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES})
  add_definitions(${HDF5_DEFINITIONS})
endif()

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
  if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

# TBB
find_package(TBB REQUIRED)
add_definitions(${TBB_DEFINITIONS})
include_directories(${TBB_INCLUDE_DIRS})
link_libraries(${TBB_LIBRARIES})

# pthreads
find_package (Threads)
link_libraries (${CMAKE_THREAD_LIBS_INIT})

# VTK
find_package(VTK)
if (VTK_FOUND)
  include_directories(${VTK_INCLUDE_DIRS})
  link_libraries(${VTK_LIBRARIES})
  add_definitions("-DUSE_VTK")
  set(mpm_vtk ${mpm_SOURCE_DIR}/src/vtk_writer.cc)
endif()

# METIS
find_package(METIS)
if (METIS_FOUND)
    include_directories(SYSTEM ${METIS_INCLUDE_DIR})
    link_libraries(${METIS_LIBRARY})
endif (METIS_FOUND)

# ParMETIS
find_package(ParMETIS)
if (PARMETIS_FOUND)
  add_definitions("-DUSE_PARMETIS")
  include_directories(${PARMETIS_INCLUDE_DIRS})
  link_libraries(${PARMETIS_LIBRARIES})
endif()


# Include directories
include_directories(BEFORE
  ${mpm_SOURCE_DIR}/include/
  ${mpm_SOURCE_DIR}/include/elements/
  ${mpm_SOURCE_DIR}/include/elements/2d
  ${mpm_SOURCE_DIR}/include/elements/3d
  ${mpm_SOURCE_DIR}/external/
  ${mpm_SOURCE_DIR}/tests/include/
)

# mpm executable
SET(mpm_src
  ${mpm_SOURCE_DIR}/src/affine_transform.cc
  ${mpm_SOURCE_DIR}/src/cell.cc
  ${mpm_SOURCE_DIR}/src/geometry.cc
  ${mpm_SOURCE_DIR}/src/hdf5_particle.cc  
  ${mpm_SOURCE_DIR}/src/io.cc
  ${mpm_SOURCE_DIR}/src/logger.cc
  ${mpm_SOURCE_DIR}/src/material.cc
  ${mpm_SOURCE_DIR}/src/mpm.cc
  ${mpm_SOURCE_DIR}/src/node.cc
  ${mpm_SOURCE_DIR}/src/particle.cc
  ${mpm_SOURCE_DIR}/src/quadrature.cc
  ${mpm_SOURCE_DIR}/src/read_mesh.cc
  ${mpm_SOURCE_DIR}/src/element.cc
)

if(MPM_BUILD_LIB)
add_library(lmpm SHARED ${mpm_src} ${mpm_vtk})
add_executable(mpm ${mpm_SOURCE_DIR}/src/main.cc)
target_link_libraries(mpm lmpm)
else()
add_executable(mpm ${mpm_SOURCE_DIR}/src/main.cc ${mpm_src} ${mpm_vtk})
endif()

# Unit test
if(MPM_BUILD_TESTING)
  SET(test_src
    ${mpm_SOURCE_DIR}/tests/test_main.cc
    ${mpm_SOURCE_DIR}/tests/cell_container_test.cc
    ${mpm_SOURCE_DIR}/tests/cell_test.cc
    ${mpm_SOURCE_DIR}/tests/factory_test.cc
    ${mpm_SOURCE_DIR}/tests/geometry_test.cc
    ${mpm_SOURCE_DIR}/tests/elements/hexahedron_element_test.cc
    ${mpm_SOURCE_DIR}/tests/elements/hexahedron_quadrature_test.cc
    ${mpm_SOURCE_DIR}/tests/elements/quadrilateral_element_test.cc
    ${mpm_SOURCE_DIR}/tests/elements/quadrilateral_gimp_element_test.cc
    ${mpm_SOURCE_DIR}/tests/elements/quadrilateral_quadrature_test.cc
    ${mpm_SOURCE_DIR}/tests/elements/triangle_element_test.cc
    ${mpm_SOURCE_DIR}/tests/elements/triangle_quadrature_test.cc
    ${mpm_SOURCE_DIR}/tests/graph_test.cc
    ${mpm_SOURCE_DIR}/tests/io_test.cc
    ${mpm_SOURCE_DIR}/tests/material/bingham_test.cc
    ${mpm_SOURCE_DIR}/tests/material/linear_elastic_test.cc
    ${mpm_SOURCE_DIR}/tests/material/newtonian_test.cc
    ${mpm_SOURCE_DIR}/tests/mesh_cell_neighbours_test.cc
    ${mpm_SOURCE_DIR}/tests/mesh_test.cc
    ${mpm_SOURCE_DIR}/tests/mpi_particle_test.cc
    ${mpm_SOURCE_DIR}/tests/mpi_transfer_particle_test.cc
    ${mpm_SOURCE_DIR}/tests/mpi_wrapper_test.cc
    ${mpm_SOURCE_DIR}/tests/mpm_explicit_usf_test.cc
    ${mpm_SOURCE_DIR}/tests/mpm_explicit_usf_unitcell_test.cc
    ${mpm_SOURCE_DIR}/tests/mpm_explicit_usl_test.cc
    ${mpm_SOURCE_DIR}/tests/mpm_explicit_usl_unitcell_test.cc
    ${mpm_SOURCE_DIR}/tests/node_container_test.cc
    ${mpm_SOURCE_DIR}/tests/node_map_test.cc
    ${mpm_SOURCE_DIR}/tests/node_test.cc
    ${mpm_SOURCE_DIR}/tests/particle_cell_crossing_test.cc
    ${mpm_SOURCE_DIR}/tests/particle_container_test.cc
    ${mpm_SOURCE_DIR}/tests/particle_test.cc
    ${mpm_SOURCE_DIR}/tests/point_in_cell_test.cc
    ${mpm_SOURCE_DIR}/tests/read_mesh_ascii_test.cc
    ${mpm_SOURCE_DIR}/tests/write_mesh_particles.cc
    ${mpm_SOURCE_DIR}/tests/write_mesh_particles_unitcell.cc
    ${mpm_SOURCE_DIR}/tests/elements/hexahedron_gimp_element_test.cc
  )
  if(MPM_BUILD_LIB)
    add_executable(mpmtest ${test_src})
    target_link_libraries(mpmtest lmpm)
    add_test(NAME mpmtest COMMAND $<TARGET_FILE:mpmtest>)
    enable_testing()
  else()
    add_executable(mpmtest ${test_src} ${mpm_src} ${mpm_vtk})
    add_test(NAME mpmtest COMMAND $<TARGET_FILE:mpmtest>)
    enable_testing()
  endif()
endif()

# Coverage
find_package(codecov)
if(ENABLE_COVERAGE)
  add_executable(mpmtest_coverage ${mpm_vtk} ${mpm_src} ${test_src})
  add_coverage(mpmtest_coverage)
endif()
