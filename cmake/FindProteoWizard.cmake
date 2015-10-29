# - Try to find Proteowizard packages
# Once done this will define
#
# This module defines
# PROTEOWIZARD_INCLUDE_DIR, where to find xercesc.h, etc.
# PROTEOWIZARD_LIBRARIES, the libraries to link against to use xercesc.
# PROTEOWIZARD_FOUND, If false, don't try to use xercesc.

function(_pwiz_find_library _name)
    find_library(${_name} NAMES ${ARGN})
    mark_as_advanced(${_name})
endfunction()

FIND_PATH(PWIZ_INCLUDE_DIR pwiz )

MESSAGE( STATUS "Found pwiz include directory: ${PWIZ_INCLUDE_DIR}" )

_pwiz_find_library(PWIZ_DATA_MSDATA_LIBRARY pwiz_data_msdata)
_pwiz_find_library(PWIZ_DATA_COMMON_LIBRARY pwiz_data_common)
_pwiz_find_library(PWIZ_UTILITY_MISC_LIBRARY pwiz_utility_misc)
_pwiz_find_library(PWIZ_UTILITY_MINIXML_LIBRARY pwiz_utility_minimxml)
_pwiz_find_library(PWIZ_DATA_MSDATA_VERSION_LIBRARY pwiz_data_msdata_version)
_pwiz_find_library(PWIZ_SHA1_LIBRARY SHA1)


set(PWIZ_LIBRARIES "")
list(APPEND PWIZ_LIBRARIES ${PWIZ_DATA_MSDATA_LIBRARY} ${PWIZ_DATA_COMMON_LIBRARY} ${PWIZ_UTILITY_MISC_LIBRARY} ${PWIZ_UTILITY_MINIXML_LIBRARY} ${PWIZ_DATA_MSDATA_VERSION_LIBRARY} ${PWIZ_SHA1_LIBRARY})

# if the include and the library are found then we have it
if (PWIZ_INCLUDE_DIR AND PWIZ_LIBRARIES)
  SET( PWIZ_FOUND 1 )
endif()

MARK_AS_ADVANCED(
  PWIZ_INCLUDE_DIR
  PWIZ_LIBRARIES
)
