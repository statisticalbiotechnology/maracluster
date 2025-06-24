# - Try to find Proteowizard packages
# Once done this will define
#
# This module defines
# PROTEOWIZARD_INCLUDE_DIR, where to find xercesc.h, etc.
# PROTEOWIZARD_LIBRARIES, the libraries to link against to use xercesc.
# PROTEOWIZARD_FOUND, If false, don't try to use xercesc.

set(PWIZ_LIBRARIES "")

macro(_pwiz_find_library _name)
  find_library(${_name} NAMES ${ARGN})
  mark_as_advanced(${_name})
  if(NOT ${_name})
    message(FATAL_ERROR "${_name} library not found")
  endif()
  list(APPEND PWIZ_LIBRARIES ${${_name}})
endmacro()

FIND_PATH(PWIZ_INCLUDE_DIR pwiz )

MESSAGE( STATUS "Found pwiz include directory: ${PWIZ_INCLUDE_DIR}" )

_pwiz_find_library(PWIZ_DATA_MSDATA_LIBRARY pwiz_data_msdata)
_pwiz_find_library(PWIZ_DATA_MSDATA_MZMLB_LIBRARY pwiz_data_msdata_mzmlb)
_pwiz_find_library(PWIZ_DATA_MSDATA_MZ5_LIBRARY pwiz_data_msdata_mz5)
_pwiz_find_library(HDF5PP_LIBRARY hdf5pp)
_pwiz_find_library(HDF5_LIBRARY hdf5)
_pwiz_find_library(PWIZ_DATA_MSDATA_CORE_LIBRARY pwiz_data_msdata_core)
_pwiz_find_library(PWIZ_DATA_COMMON_LIBRARY pwiz_data_common)
_pwiz_find_library(PWIZ_UTILITY_MISC_LIBRARY pwiz_utility_misc)
_pwiz_find_library(PWIZ_UTILITY_MINIXML_LIBRARY pwiz_utility_minimxml)
_pwiz_find_library(PWIZ_DATA_MSDATA_VERSION_LIBRARY pwiz_data_msdata_version)
_pwiz_find_library(PWIZ_SHA1_LIBRARY SHA1)

#/home/matthew/build-32bit/tools/proteowizard/lib/baf2sql_c.dll
#/home/matthew/build-32bit/tools/proteowizard/lib/cdt.dll
#/home/matthew/build-32bit/tools/proteowizard/lib/MassLynxRaw.dll
if(MSVC)
  _pwiz_find_library(SQLITE3PP_LIBRARY sqlite3pp)
  _pwiz_find_library(SQLITE3_LIBRARY sqlite3)
  if (VENDOR_SUPPORT)
    _pwiz_find_library(PWIZ_DATA_VENDOR_LIBRARY pwiz_data_vendor_readers)
    _pwiz_find_library(PWIZ_DATA_VENDOR_THERMO_LIBRARY pwiz_reader_thermo)
    _pwiz_find_library(PWIZ_DATA_VENDOR_THERMO_API_LIBRARY pwiz_vendor_api_thermo)
    _pwiz_find_library(PWIZ_DATA_VENDOR_SHIMADZU_LIBRARY pwiz_reader_shimadzu)
    _pwiz_find_library(PWIZ_DATA_VENDOR_SHIMADZU_API_LIBRARY pwiz_vendor_api_shimadzu)
    _pwiz_find_library(PWIZ_DATA_VENDOR_UIMF_LIBRARY pwiz_reader_uimf)
    _pwiz_find_library(PWIZ_DATA_VENDOR_UIMF_API_LIBRARY pwiz_vendor_api_uimf)
    _pwiz_find_library(PWIZ_DATA_VENDOR_UNIFI_LIBRARY pwiz_reader_unifi)
    _pwiz_find_library(PWIZ_DATA_VENDOR_UNIFI_API_LIBRARY pwiz_vendor_api_unifi)
    _pwiz_find_library(PWIZ_DATA_VENDOR_AGILENT_LIBRARY pwiz_reader_agilent)
    _pwiz_find_library(PWIZ_DATA_VENDOR_AGILENT_API_LIBRARY pwiz_vendor_api_agilent)
    _pwiz_find_library(PWIZ_DATA_VENDOR_WATERS_LIBRARY pwiz_reader_waters)
    _pwiz_find_library(PWIZ_DATA_VENDOR_WATERS_API_LIBRARY MassLynxRaw)
    _pwiz_find_library(PWIZ_DATA_VENDOR_WATERS_IMS_API_LIBRARY cdt)
    _pwiz_find_library(PWIZ_DATA_VENDOR_BRUKER_LIBRARY pwiz_reader_bruker)
    _pwiz_find_library(PWIZ_DATA_VENDOR_BRUKER_BAF2SQL_API_LIBRARY baf2sql_c)
    _pwiz_find_library(PWIZ_DATA_VENDOR_BRUKER_TIMSDATA_API_LIBRARY timsdata)
    _pwiz_find_library(PWIZ_DATA_VENDOR_BRUKER_API_LIBRARY pwiz_vendor_api_bruker)
    _pwiz_find_library(PWIZ_DATA_VENDOR_ABI_T2D_LIBRARY pwiz_reader_abi_t2d)
    _pwiz_find_library(PWIZ_DATA_VENDOR_ABI_LIBRARY pwiz_reader_abi)
    _pwiz_find_library(PWIZ_DATA_VENDOR_ABI_API_LIBRARY pwiz_vendor_api_abi)
    if (NOT CMAKE_GENERATOR MATCHES ".*Win64.*")
      _pwiz_find_library(PWIZ_DATA_VENDOR_ABI_T2D_API_LIBRARY pwiz_vendor_api_abi_t2d)
    endif()
  endif (VENDOR_SUPPORT)
endif(MSVC)

# if the include and the library are found then we have it
if (PWIZ_INCLUDE_DIR AND PWIZ_LIBRARIES)
  SET( PWIZ_FOUND 1 )
endif()

MARK_AS_ADVANCED(
  PWIZ_INCLUDE_DIR
  PWIZ_LIBRARIES
)
