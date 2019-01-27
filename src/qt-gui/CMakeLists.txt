cmake_minimum_required(VERSION 2.8)
project(MARACLUSTER_GUI)
set(MARACLUSTER_GUI_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

if (UNIX)
  set(CMAKE_INSTALL_RPATH "$ORIGIN")
  set(MARACLUSTER_BIN_DESTINATION lib/maracluster)
endif (UNIX)

if (APPLE)
  set(MARACLUSTER_BIN_DESTINATION maracluster-gui.app/Contents/MacOS)
endif (APPLE)

# get the maracluster executable
set(SKIP_MAIN_CPACK ON)
set(MARACLUSTER_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../..")
add_subdirectory(${MARACLUSTER_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/main)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/main/src)

# variables set by add_subdirectory are not transferred to the parent scope
set(CMAKE_MODULE_PATH ${MARACLUSTER_SOURCE_DIR}/cmake)
include("${MARACLUSTER_SOURCE_DIR}/CommonCMake.txt")

set(CMAKE_INSTALL_DO_STRIP ON)

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

set(Qt_libs Qt5::Core Qt5::Gui Qt5::Widgets)

if (APPLE)
  find_package(Qt5DBus CONFIG REQUIRED)
  find_package(Qt5PrintSupport CONFIG REQUIRED)
  
  set(ICON_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/icon/maracluster.icns)
  set_source_files_properties(${ICON_CONFIG} PROPERTIES
         MACOSX_PACKAGE_LOCATION "Resources")
  list(APPEND Qt_libs Qt5::DBus Qt5::PrintSupport)
elseif (UNIX)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/icon/hicolor DESTINATION share/icons)
elseif (MSVC)
  set(ICON_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/packaging/maracluster-gui.rc)
endif (APPLE)

# Tell CMake to create the executable
add_executable(maracluster-gui MACOSX_BUNDLE WIN32 ${ICON_CONFIG} ${maraclustergui_SRCS})

target_link_libraries(maracluster-gui ${Qt_libs})

install(TARGETS maracluster-gui EXPORT MARACLUSTER BUNDLE DESTINATION . RUNTIME DESTINATION lib/maracluster) # Important to use relative path here (used by CPack)!

### INSTALL DEPENDENCIES (PLATFORM DEPENDENT) ###

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
elseif (UNIX)
  set(CMAKE_FIND_LIBRARY_SUFFIXES_BACKUP ${CMAKE_FIND_LIBRARY_SUFFIXES})
  set(CMAKE_FIND_LIBRARY_SUFFIXES .so) # only look for shared libraries!
  
  # these two packages (DBus and XcbQpa) are required by the Xcb Plugin
  find_package(Qt5DBus CONFIG REQUIRED)
  find_library(QtXcbQpa_LIBRARY_FULL NAMES libQt5XcbQpa.so )
  
  # the c++ standard, freetype and ICU libraries are often not compatible with older versions
  find_package(ICU REQUIRED COMPONENTS uc i18n data)
  find_package(Freetype REQUIRED)
  find_library(STDCPP_LIBRARY NAMES libstdc++.so.6 )
  
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_BACKUP})  
  
  macro(_qt5_get_shared_library lib)
    get_target_property(Qt${lib}_location Qt5::${lib} LOCATION)
    get_target_property(Qt${lib}_soname Qt5::${lib} IMPORTED_SONAME_RELEASE)
    get_filename_component(Qt${lib}_path ${Qt${lib}_location} PATH)
    set(Qt${lib}_so_path "${Qt${lib}_path}/${Qt${lib}_soname}")
    set(QT_SO_LIBRARIES ${QT_SO_LIBRARIES} ${Qt${lib}_so_path} ${Qt${lib}_location})
  endmacro()
  
  _qt5_get_shared_library(Core)
  _qt5_get_shared_library(Widgets)
  _qt5_get_shared_library(Gui)
  _qt5_get_shared_library(DBus)
  
  foreach (LIBRARY ${ICU_LIBRARIES} ${FREETYPE_LIBRARIES} ${QtXcbQpa_LIBRARY_FULL} ${STDCPP_LIBRARY})
    get_filename_component(LIBRARY_RESOLVED ${LIBRARY} REALPATH)
    list(APPEND LIBRARIES_RESOLVED ${LIBRARY_RESOLVED})
    string(REGEX REPLACE "(.*\.so.[0-9]+).*" "\\1" LIBRARY_MAIN_VERSION "${LIBRARY_RESOLVED}")
    list(APPEND LIBRARIES_RESOLVED ${LIBRARY_MAIN_VERSION})
  endforeach (LIBRARY)
  
  message( STATUS "Installing shared libraries: ${LIBRARIES_RESOLVED};${QT_SO_LIBRARIES}")
  
  install(PROGRAMS ${QT_SO_LIBRARIES} ${LIBRARIES_RESOLVED} DESTINATION lib/maracluster)
  
  # the libqxcb.so shared library has to be located in a folder named "platforms" relative to the executable
  get_target_property(QtXcbPlugin_location Qt5::QXcbIntegrationPlugin LOCATION)
  get_filename_component(QtXcbPlugin_filename ${QtXcbPlugin_location} NAME)
  
  file(COPY ${QtXcbPlugin_location} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  set(QtXcbPlugin_copied_location ${CMAKE_CURRENT_BINARY_DIR}/${QtXcbPlugin_filename}) # runs into permissions error
  find_program (PATCHELF_BINARY NAMES patchelf)
  message( STATUS "Found patchelf binary: ${PATCHELF_BINARY}")
  install(CODE "execute_process(COMMAND ${PATCHELF_BINARY} --remove-rpath ${QtXcbPlugin_copied_location})") # gives "stat: No such file or directory" error on earlier versions of patchelf. This can be ignored.
  install(CODE "execute_process(COMMAND ${PATCHELF_BINARY} --force-rpath --set-rpath '$ORIGIN/..' ${QtXcbPlugin_copied_location})")
  
  install(PROGRAMS ${QtXcbPlugin_copied_location} DESTINATION lib/maracluster/platforms)
  
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/packaging/maracluster.desktop.cmake ${CMAKE_CURRENT_BINARY_DIR}/packaging/maracluster.desktop)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/packaging/maracluster.desktop DESTINATION share/applications)
elseif (MSVC)
  get_target_property(QtCore_location Qt5::Core LOCATION)
  get_target_property(QtWidgets_location Qt5::Widgets LOCATION)
  get_target_property(QtGui_location Qt5::Gui LOCATION)
  
  set(DLL_LIBRARIES ${QtCore_location} ${QtGui_location} ${QtWidgets_location})
  
  message( STATUS "    Win32 DLLs ${DLL_LIBRARIES}." )  
  set( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS ${DLL_LIBRARIES} )
  include( InstallRequiredSystemLibraries )
endif (APPLE)

set(PACKAGE_NICKNAME "maracluster-gui")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Clustering spectra")

include("CPack.txt")
include(CPack)
