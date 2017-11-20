# - Try to find valgrind
# Once done this will define
#  VALGRIND_FOUND       - System has valgrind
#  VALGRIND_EXECUTABLE  - The executable location

FIND_PROGRAM(VALGRIND_EXECUTABLE NAME valgrind
             HINTS $ENV{VALGRIND_HOME}
                   "/usr/bin/"
                   "/usr/local/bin")

# handle the QUIETLY and REQUIRED arguments and set VALGRIND_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VALGRIND DEFAULT_MSG VALGRIND_EXECUTABLE)

MARK_AS_ADVANCED(VALGRIND_EXECUTABLE)
