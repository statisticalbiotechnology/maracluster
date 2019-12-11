#############################################################################
# NAME AND VERSION
###############################################################################
cmake_minimum_required(VERSION 2.8)
project(MARACLUSTER)
if (NOT DEFINED MARACLUSTER_SOURCE_DIR)
  set(MARACLUSTER_SOURCE_DIR ${CMAKE_SOURCE_DIR})
endif ()
include("CommonCMake.txt")

###############################################################################
# PREPARING TO INSTALL
###############################################################################

my_set(CMAKE_BUILD_TYPE "Debug" "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel.")
#my_set(CMAKE_PREFIX_PATH "../" "Default path to packages")

option(TARBALL_BUILD "Create a tarball instead of an installer package." OFF)
if(TARBALL_BUILD)
  add_definitions(-DTARBALL_BUILD) # triggers ifdef statements in the C++ code
endif(TARBALL_BUILD)

option(FINGERPRINT_FILTER "Use pre-filtering of spectrum pairs by their fragment fingerprints." OFF)
if(FINGERPRINT_FILTER)
  add_definitions(-DFINGERPRINT_FILTER) # triggers ifdef statements in the C++ code
endif(FINGERPRINT_FILTER)

option(DOT_PRODUCT "Use dot product similarity measure instead of p values." OFF)
if(DOT_PRODUCT)
  add_definitions(-DDOT_PRODUCT) # triggers ifdef statements in the C++ code
endif(DOT_PRODUCT)

option(SINGLE_LINKAGE "Use single instead of complete linkage for clustering." OFF)
if(SINGLE_LINKAGE)
  add_definitions(-DSINGLE_LINKAGE) # triggers ifdef statements in the C++ code
endif(SINGLE_LINKAGE)

option(VENDOR_SUPPORT "Add support for reading vendor RAW files using Proteowizard." OFF)
if(VENDOR_SUPPORT)
  add_definitions(-DVENDOR_SUPPORT) # triggers ifdef statements in the C++ code
endif(VENDOR_SUPPORT)

# PRINT VARIBALES TO STDOUT
MESSAGE( STATUS )
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS "Building MARACLUSTER:" )
MESSAGE( STATUS "change a configuration variable with: cmake -D<Variable>=<Value>" )
MESSAGE( STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}" )
MESSAGE( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}" )
MESSAGE( STATUS "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}, ${CMAKE_MODULE_PATH}" )
MESSAGE( STATUS "TARGET_ARCH = ${TARGET_ARCH}" )
MESSAGE( STATUS "TOOL CHAIN FILE = ${CMAKE_TOOLCHAIN_FILE}")
MESSAGE( STATUS "PROFILING = ${PROFILING}")
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)

MESSAGE( STATUS "MARACLUSTER options:" )
MESSAGE( STATUS "FINGERPRINT_FILTER = ${FINGERPRINT_FILTER}")
MESSAGE( STATUS "DOT_PRODUCT = ${DOT_PRODUCT}")
MESSAGE( STATUS "SINGLE_LINKAGE = ${SINGLE_LINKAGE}")
MESSAGE( STATUS "VENDOR_SUPPORT = ${VENDOR_SUPPORT}")
MESSAGE( STATUS
"-------------------------------------------------------------------------------"
)
MESSAGE( STATUS )


# Detect machine architecture, on UNIX:
get_arch("${TARGET_ARCH}")

# STORE NEWLY SET VARIABLES IN *.h.cmake FILES
load_config_files("${CMAKE_CURRENT_SOURCE_DIR}")

###############################################################################
# COMPILING
###############################################################################
add_subdirectory(src)

###############################################################################
# TESTING
###############################################################################

# Enabling system level tests (Ctest)
#enable_testing()
# Scheduling system level tests
#add_subdirectory(data/system_tests/percolator)
# Scheduling unit level tests
#if(GOOGLE_TEST)
#  add_subdirectory(data/unit_tests/percolator)
#endif()

###############################################################################
# INSTALLING
###############################################################################


if( MINGW )
  message( STATUS "    Installing system-libraries: MinGW DLLs." )
  #set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${MINGW_PATH}/bin/libxerces-c-3-1.dll)
  include( InstallRequiredSystemLibraries )
endif( MINGW )

if( MSVC )
  message( STATUS "    Installing system-libraries: Win32 DLLs." )
  set( CMAKE_INSTALL_OPENMP_LIBRARIES TRUE )
  if (VENDOR_SUPPORT)
    set( DLL_LIBRARIES "${PWIZ_DATA_VENDOR_BRUKER_BAF2SQL_API_LIBRARY};${PWIZ_DATA_VENDOR_WATERS_API_LIBRARY};${PWIZ_DATA_VENDOR_WATERS_IMS_API_LIBRARY}" )
    STRING(REPLACE ".lib" ".dll" DLL_LIBRARIES "${DLL_LIBRARIES}")
    set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${DLL_LIBRARIES} )
  endif(VENDOR_SUPPORT)
  include( InstallRequiredSystemLibraries )
endif( MSVC )


# ###############################################################################
# # PACKAGING
# ###############################################################################
# # PACKAGING OPTIONS: GENERAL

if(VENDOR_SUPPORT)
  set(NICKNAME_SUFFIX "-vendor-support")
endif(VENDOR_SUPPORT)

set(PACKAGE_NICKNAME "maracluster${NICKNAME_SUFFIX}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Clustering spectra")
set(CPACK_RPM_PACKAGE_DEPENDS "libc6, libgcc1")

if (NOT SKIP_MAIN_CPACK)
  include("CPack.txt")
  include(CPack)
else (NOT SKIP_MAIN_CPACK)
  message(STATUS "Skipping main CPack")
endif (NOT SKIP_MAIN_CPACK)
