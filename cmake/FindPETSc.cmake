#################################################################
# Try to find PETSc                                             #
#                                                               #
# Once done this will define:                                   #
#  PETSC_FOUND         - system has PETSc                       #
#  PETSC_DIR           - PETSc directory                        #
#  PETSC_ARCH          - PETSc architecture                     #
#  PETSC_INC           - PETSc include directory                #
#  PETSC_LIB           - PETSc library (static or dynamic)      #
#  PETSC_VARIABLES     - Content of PETSc 'petscvariables' file #
#                                                               #
#  PETSC_MUMPS         - Was PETSc compiled with MUMPS?         #
#  PETSC_MUMPS_INC     - PETSc MUMPS include file               #
#  PETSC_MUMPS_LIB     - PETSc MUMPS libraries                  #
#                                                               #
#  PETSC_SCALAPACK     - Was PETSc compiled with ScaLAPACK?     #
#  PETSC_SCALAPACK_LIB - PETSc ScaLAPACK libraries              #
#                                                               #
#  PETSC_PARMETIS      - Was PETSc compiled with ParMETIS?      #
#  PETSC_PARMETIS_LIB  - PETSc ParMETIS libraries               #
#                                                               #
#  PETSC_METIS         - Was PETSc compiled with Metis?         #
#  PETSC_METIS_LIB     - PETSc METIS libraries                  #
#                                                               #
#  PETSC_MPIUNI        - Was PETSc compiled with MPIUNI?        #
#  PETSC_MPIUNI_INC    - MPIUNI include file                    #
#                                                               #
# Usage:                                                        #
#  find_package(PETSc)                                          #
#                                                               #
# Setting these changes the behavior of the search              #
#  PETSC_DIR           - PETSc directory                        #
#  PETSC_ARCH          - PETSc architecture                     #
#################################################################

## Try to set PETSC_DIR and PETSC_ARCH ##
#########################################
if(NOT DEFINED PETSC_DIR)
  set(PETSC_DIR $ENV{PETSC_DIR})
endif()
if(NOT DEFINED PETSC_ARCH)
  set(PETSC_ARCH $ENV{PETSC_ARCH})
endif()

## Includes ##
##############
if(EXISTS "${PETSC_DIR}/include" AND
   EXISTS "${PETSC_DIR}/${PETSC_ARCH}/include")
  set(PETSC_INC "${PETSC_DIR}/include" "${PETSC_DIR}/${PETSC_ARCH}/include")
else()
  message(SEND_ERROR "PETSc includes not found")
endif()

## Library ##
#############
if(EXISTS "${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so")
  set(PETSC_LIB "${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so")
elseif(EXISTS "${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.a")
  set(PETSC_LIB "${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.a")
else()
  message(SEND_ERROR "PETSc library not found")
endif()

## PETSc variables ##
#####################
if(EXISTS ${PETSC_DIR}/${PETSC_ARCH}/conf/petscvariables)
  file(STRINGS ${PETSC_DIR}/${PETSC_ARCH}/conf/petscvariables
    PETSC_VARIABLES NEWLINE_CONSUME)
elseif(EXISTS ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables)
  file(STRINGS ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables
    PETSC_VARIABLES NEWLINE_CONSUME)
else()
  message(SEND_ERROR "PETSc variables not found")
endif()

## MUMPS ##
###########
if(EXISTS ${PETSC_DIR}/${PETSC_ARCH}/lib/libmumps_common.a)
  set(PETSC_MUMPS TRUE)
  set(PETSC_MUMPS_INC ${PETSC_INC})
  set(PETSC_MUMPS_LIB
    ${PETSC_DIR}/${PETSC_ARCH}/lib/libcmumps.a
    ${PETSC_DIR}/${PETSC_ARCH}/lib/libdmumps.a
    ${PETSC_DIR}/${PETSC_ARCH}/lib/libsmumps.a
    ${PETSC_DIR}/${PETSC_ARCH}/lib/libzmumps.a
    ${PETSC_DIR}/${PETSC_ARCH}/lib/libmumps_common.a
    ${PETSC_DIR}/${PETSC_ARCH}/lib/libpord.a)

  if(EXISTS ${PETSC_DIR}/${PETSC_ARCH}/lib/libmpiseq.a)
    set(PETSC_MUMPS_SEQ ${PETSC_DIR}/${PETSC_ARCH}/lib/libmpiseq.a)
  else()
    set(PETSC_MUMPS_SEQ "")
  endif()

else()
  set(PETSC_MUMPS FALSE)
  set(PETSC_MUMPS_INC "")
  set(PETSC_MUMPS_LIB "")
  set(PETSC_MUMPS_SEQ "")
endif()

## ScaLAPACK ##
###############
if(EXISTS ${PETSC_DIR}/${PETSC_ARCH}/lib/libscalapack.a)
  set(PETSC_SCALAPACK TRUE)
  set(PETSC_SCALAPACK_LIB ${PETSC_DIR}/${PETSC_ARCH}/lib/libscalapack.a)
else()
  set(PETSC_SCALAPACK FALSE)
  set(PETSC_SCALAPACK_LIB "")
endif()

## ParMETIS ##
##############
if(EXISTS ${PETSC_DIR}/${PETSC_ARCH}/lib/libparmetis.so)
  set(PETSC_PARMETIS TRUE)
  set(PETSC_PARMETIS_LIB ${PETSC_DIR}/${PETSC_ARCH}/lib/libparmetis.so)
else()
  set(PETSC_PARMETIS FALSE)
  set(PETSC_PARMETIS_LIB "")
endif()

## METIS ##
############
if(EXISTS ${PETSC_DIR}/${PETSC_ARCH}/lib/libmetis.so)
  set(PETSC_METIS TRUE)
  set(PETSC_METIS_LIB ${PETSC_DIR}/${PETSC_ARCH}/lib/libmetis.so)
else()
  set(PETSC_METIS FALSE)
  set(PETSC_METIS_LIB "")
endif()

## MPI Uni ##
#############
string(REGEX MATCH "MPI_IS_MPIUNI = [^\n\r]*" PETSC_MPIUNI ${PETSC_VARIABLES})
if(PETSC_MPIUNI)
  string(REPLACE "MPI_IS_MPIUNI = " "" PETSC_MPIUNI ${PETSC_MPIUNI})
  string(STRIP ${PETSC_MPIUNI} PETSC_MPIUNI)
  string(COMPARE EQUAL ${PETSC_MPIUNI} "1" PETSC_MPIUNI)

  string(REGEX MATCH "MPI_INCLUDE = -I[^\n\r]*"
    PETSC_MPIUNI_INC ${PETSC_VARIABLES})
  string(REPLACE "MPI_INCLUDE = -I" "" PETSC_MPIUNI_INC ${PETSC_MPIUNI_INC})
  string(STRIP ${PETSC_MPIUNI_INC} PETSC_MPIUNI_INC)
else()
  set(PETSC_MPIUNI     FALSE)
  set(PETSC_MPIUNI_INC "")
endif()

## CMake check and done ##
##########################
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc
  "PETSc could not be found: be sure to set PETSC_DIR and PETSC_ARCH in your environment variables"
  PETSC_LIB PETSC_INC PETSC_DIR PETSC_ARCH)
