#
# Findgalaxy
# ----------# 
#
find_path(GXY_INCLUDE_DIR NAMES galaxy.h PATH_SUFFIXES gxy)
find_library(GXY_LIBRARIES NAMES gxy_framework gxy_data)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GXY DEFAULT_MSG GXY_INCLUDE_DIR GXY_LIBRARIES)

if(DEFINED GXY_LIBRARY_PATH)
  find_library(GXY_LIBRARY_FILE_TEST gxy_data libgxy_framework libgxy_module_data libgxy_data libgxy_renderer HINTS ${DEFAULT_LIBRARY_PATH})
  message(STATUS "Undefined case (blank string): Searching for librarary files for gxy_data in ${DEFAULT_LIBRARY_PATH}")
endif(DEFINED GXY_LIBRARY_PATH)

set(GXY_LIBRARY_DIR ${GXY_ROOT_DIR}/lib64)
set(GALAXY_LIBRARIES gxy_framework gxy_renderer gxy_data gxy_ospray)
#set(GXY_RENDERING_LIBRARIES libgxy_renderer.so)

mark_as_advanced(GXY_INCLUDE_DIR)
mark_as_advanced(GXY_LIBRARIES)
