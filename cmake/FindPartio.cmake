# - Try to find Partio
# Once done this will define
#
#  PARTIO_FOUND        - system has Partio
#  PARTIO_INCLUDE_DIRS - include directories for Partio
#  PARTIO_LIBRARIES    - libraries for Partio
#=============================================================================


message(STATUS "Checking for package 'Partio'")

# Check for header file
find_path(PARTIO_INCLUDE_DIRS Partio.h
  HINTS ${PARTIO_ROOT}/src/lib
  PATH_SUFFIXES partio
  DOC "Directory where the Partio header files are located"
  )

find_library(PARTIO_LIBRARY partio
  HINTS ${PARTIO_ROOT}/src/lib
  NO_DEFAULT_PATH
  DOC "Directory where the Partio library is located"
  )

find_library(PARTIO_LIBRARY partio
  DOC "Directory where the Partio library is located"
  )

set(PARTIO_LIBRARIES ${PARTIO_LIBRARY})


find_package_handle_standard_args(Partio
                                  "Partio could not be found/configured."
                                  PARTIO_LIBRARIES
                                  PARTIO_INCLUDE_DIRS
                                )
