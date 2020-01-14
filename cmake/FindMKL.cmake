#
# FindMKL++
# ----------
#
# Try to find the Intel Math Kernel Library

if(NOT MKL_FOUND)
  set(MKL_INCLUDE_PATH "$ENV{MKLROOT}/include")
  set(MKL_LIBRARY_PATH "$ENV{MKLROOT}/lib/intel64")

  find_path(MKL_INCLUDE_DIR mkl.h HINTS ${MKL_INCLUDE_PATH} ENV CPATH)
  find_library(MKL_CORE mkl_core HINTS ${MKL_LIBRARY_PATH} ENV LD_LIBRARY_PATH)
  find_library(MKL_THREAD mkl_intel_thread HINTS ${MKL_LIBRARY_PATH} ENV LD_LIBRARY_PATH)
  find_library(MKL_LP mkl_intel_lp64 HINTS ${MKL_LIBRARY_PATH} ENV LD_LIBRARY_PATH)

  set(MKL_LIBRARIES ${MKL_LP} ${MKL_THREAD} ${MKL_CORE})

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDE_DIR MKL_LIBRARIES)

  mark_as_advanced(MKL_INCLUDE_DIR MKL_LIBRARIES)
endif()
