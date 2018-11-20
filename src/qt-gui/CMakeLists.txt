cmake_minimum_required(VERSION 3.5.0)

project(maracluster-gui)
set(MARACLUSTER_SOURCE_DIR ${CMAKE_SOURCE_DIR})

set(CMAKE_MODULE_PATH ${MARACLUSTER_SOURCE_DIR}/cmake)

set(CMAKE_INSTALL_DO_STRIP ON)

set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "05")
set(CPACK_PACKAGE_VERSION_PATCH "0")

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)
# Create code from a list of Qt designer ui files
set(CMAKE_AUTOUIC ON)

find_package(Qt5Core CONFIG REQUIRED)
find_package(Qt5Gui CONFIG REQUIRED)
find_package(Qt5Widgets CONFIG REQUIRED)

# Populate a CMake variable with the sources
set(maraclustergui_SRCS
    mainwidget.cpp
    main.cpp
)

if (APPLE)
  set(CMAKE_MACOSX_RPATH 1)
  find_package(Qt5DBus CONFIG REQUIRED)
  find_package(Qt5PrintSupport CONFIG REQUIRED)
  
  set(ICON_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/icon/maracluster.icns)
  set_source_files_properties(${ICON_CONFIG} PROPERTIES
         MACOSX_PACKAGE_LOCATION "Resources")
  set(Qt_libs Qt5::Core Qt5::Gui Qt5::Widgets Qt5::DBus Qt5::PrintSupport)
elseif (MSVC)
  set(ICON_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/packaging/maracluster-gui.rc)
  set(Qt_libs Qt5::Core Qt5::Gui Qt5::Widgets)
elseif (UNIX)
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib/maracluster")
  set(Qt_libs Qt5::Core Qt5::Gui Qt5::Widgets)
endif (APPLE)

# Tell CMake to create the executable
add_executable(maracluster-gui MACOSX_BUNDLE WIN32 ${ICON_CONFIG} ${maraclustergui_SRCS})

target_link_libraries(maracluster-gui ${Qt_libs})

install(TARGETS maracluster-gui EXPORT MARACLUSTER BUNDLE DESTINATION . RUNTIME DESTINATION bin) # Important to use relative path here (used by CPack)!


### INSTALL DEPENDENCIES (PLATFORM DEPENDENT) ###

macro(_qt5_get_shared_library lib)
  get_target_property(Qt${lib}_location Qt5::${lib} LOCATION)
  get_target_property(Qt${lib}_soname Qt5::${lib} IMPORTED_SONAME_RELEASE)
  get_filename_component(Qt${lib}_path ${Qt${lib}_location} PATH)
  set(Qt${lib}_so_path "${Qt${lib}_path}/${Qt${lib}_soname}")
  set(QT_SO_LIBRARIES ${QT_SO_LIBRARIES} ${Qt${lib}_so_path} ${Qt${lib}_location})
endmacro()

if (UNIX AND NOT APPLE)
  find_package(ICU REQUIRED COMPONENTS uc i18n data)
  
  find_package(Freetype REQUIRED)
  find_library (
    FREETYPE_LIBRARIES 
    NAMES libfreetype.so.6
    )
  
  message( STATUS "${ICU_VERSION} ${FREETYPE_VERSION_STRING}")
  
  _qt5_get_shared_library(Core)
  _qt5_get_shared_library(Widgets)
  _qt5_get_shared_library(Gui)
  
  # these two packages (DBus and XcbQpa) are required by the Xcb Plugin
  find_package(Qt5DBus CONFIG REQUIRED)
  _qt5_get_shared_library(DBus)
  find_library (
    QtXcbQpa_LIBRARY_FULL
    NAMES libQt5XcbQpa.so.${Qt5Core_VERSION_STRING}
    )
  find_library (
    QtXcbQpa_LIBRARY
    NAMES libQt5XcbQpa.so.5
    )
  set(QT_SO_LIBRARIES ${QT_SO_LIBRARIES} ${QtXcbQpa_LIBRARY_FULL} ${QtXcbQpa_LIBRARY})
  
  message( STATUS "Installing shared libraries: ${ICU_LIBRARIES};${FREETYPE_LIBRARIES};${QT_SO_LIBRARIES}")
  
  install(PROGRAMS ${QT_SO_LIBRARIES} ${FREETYPE_LIBRARIES} ${ICU_LIBRARIES} DESTINATION lib/maracluster)
  
  get_target_property(QtXcbPlugin_location Qt5::QXcbIntegrationPlugin LOCATION)
  
  get_filename_component(QtXcbPlugin_filename ${QtXcbPlugin_location} NAME)
  install(PROGRAMS ${QtXcbPlugin_location} DESTINATION bin/platforms)
  
  install(CODE "execute_process(COMMAND patchelf --set-rpath '$ORIGIN/../../lib/maracluster' ${CMAKE_INSTALL_PREFIX}/bin/platforms/${QtXcbPlugin_filename})")
  
  # Problem with static linking of Qt: https://bugreports.qt.io/browse/QTBUG-38913
  # obtained dependencies from ${_qt5Core_install_prefix}/lib/pkfconfig/Qt5Core.pc
  #set(Qt5Core_lib_deps pthread m z icui18n icuuc icudata ${_qt5Core_install_prefix}/lib/libqtpcre2.a dl)
  #set(Qt5Gui_lib_deps m z icui18n icuuc icudata ${_qt5Core_install_prefix}/lib/libqtpcre2.a dl pthread GL png12 ${_qt5Core_install_prefix}/lib/libqtharfbuzz.a z)
  #set(Qt5Widgets_lib_deps GL png12 ${_qt5Core_install_prefix}/lib/libqtharfbuzz.a m z icui18n icuuc icudata ${_qt5Core_install_prefix}/lib/libqtpcre2.a dl pthread GL)
endif (UNIX AND NOT APPLE)

if (MSVC)
  get_target_property(QtCore_location Qt5::Core LOCATION)
  get_target_property(QtWidgets_location Qt5::Widgets LOCATION)
  get_target_property(QtGui_location Qt5::Gui LOCATION)
  
  set(DLL_LIBRARIES ${QtCore_location} ${QtGui_location} ${QtWidgets_location})
  
  message( STATUS "    Win32 DLLs ${DLL_LIBRARIES}." )  
  set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${DLL_LIBRARIES} )
  include( InstallRequiredSystemLibraries )
endif (MSVC)

if (APPLE)
  get_target_property(QtCore_location Qt5::Core LOCATION)
  get_filename_component(QtCore_path ${QtCore_location} PATH)
  
  get_target_property(QtWidgets_location Qt5::Widgets LOCATION)
  get_filename_component(QtWidgets_path ${QtWidgets_location} PATH)
  
  get_target_property(QtGui_location Qt5::Gui LOCATION)
  get_filename_component(QtGui_path ${QtGui_location} PATH)
  
  get_target_property(QtDBus_location Qt5::DBus LOCATION)
  get_filename_component(QtDBus_path ${QtDBus_location} PATH)

  get_target_property(QtPrintSupport_location Qt5::PrintSupport LOCATION)
  get_filename_component(QtPrintSupport_path ${QtPrintSupport_location} PATH)
  
  get_target_property(QtCocoaPlugin_location Qt5::QCocoaIntegrationPlugin LOCATION)
  get_filename_component(QtCocoaPlugin_name ${QtCocoaPlugin_location} NAME)
  install(PROGRAMS ${QtCocoaPlugin_location} DESTINATION maracluster-gui.app/Contents/plugins/platforms)
  set(QtCocoaPlugin_lib \${CMAKE_INSTALL_PREFIX}/maracluster-gui.app/Contents/plugins/platforms/${QtCocoaPlugin_name})
  
  get_target_property(QtMacStylePlugin_location Qt5::QMacStylePlugin LOCATION)
  get_filename_component(QtMacStylePlugin_name ${QtMacStylePlugin_location} NAME)
  install(PROGRAMS ${QtMacStylePlugin_location} DESTINATION maracluster-gui.app/Contents/plugins/styles)
  set(QtMacStylePlugin_lib \${CMAKE_INSTALL_PREFIX}/maracluster-gui.app/Contents/plugins/styles/${QtMacStylePlugin_name})
  
  install(CODE "
    include(BundleUtilities)
    fixup_bundle(\${CMAKE_INSTALL_PREFIX}/maracluster-gui.app \"${QtCocoaPlugin_lib};${QtMacStylePlugin_lib}\" \"${QtCore_path};${QtGui_path};${QtWidgets_path};${QtDBus_path};${QtPrintSupport_path}\")
  ")
endif (APPLE)

set(PACKAGE_NICKNAME "maracluster-gui")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Clustering spectra")

include("CPack.txt")
include(CPack)
