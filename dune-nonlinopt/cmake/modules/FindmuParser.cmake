# - Try to find muParser
# Once done this will define
#  MUPARSER_FOUND - System has muParser
#  MUPARSER_INCLUDE_DIRS - The muParser include directories
#  MUPARSER_LIBRARIES - The libraries needed to use muParser
#  MUPARSER_DEFINITIONS - Compiler switches required for using muParser

find_package(PkgConfig)
pkg_check_modules(PC_MUPARSER QUIET muparser)
set(MUPARSER_DEFINITIONS ${PC_MUPARSER_CFLAGS_OTHER})

find_path(MUPARSER_INCLUDE_DIR muParser.h
          HINTS ${PC_MUPARSER_INCLUDEDIR} ${PC_MUPARSER_INCLUDE_DIRS})

find_library(MUPARSER_LIBRARY NAMES muparser libmuparser
             HINTS ${PC_MUPARSER_LIBDIR} ${PC_MUPARSER_LIBRARY_DIRS})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set muParser_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(muParser DEFAULT_MSG MUPARSER_LIBRARY MUPARSER_INCLUDE_DIR)

mark_as_advanced(MUPARSER_INCLUDE_DIR MUPARSER_LIBRARY)

set(MUPARSER_LIBRARIES ${MUPARSER_LIBRARY})
set(MUPARSER_INCLUDE_DIRS ${MUPARSER_INCLUDE_DIR})
