:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::                                                                           :::
::: tested on a Windows 10 64-bit installation with VS2015 Community          :::
:::                                                                           :::
::: Support for Thermo RAW Files seems to be limited to 64-bit builds.        :::
:::                                                                           :::
::: To support reading Waters RAW Files, Visual C++ runtime is necessary to   :::
::: run http://www.microsoft.com/en-us/download/details.aspx?id=30679         :::
:::                                                                           :::
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@echo off

set SRC_DIR=%~dp0..\..\..\
set BUILD_DIR=%SRC_DIR%\build\win32
set RELEASE_DIR=%SRC_DIR%\release\win32
set BUILD_TYPE=Release
set NO_GUI="false"

:parse
IF "%~1"=="" GOTO endparse
IF "%~1"=="-s" (set SRC_DIR=%~2)
IF "%~1"=="-b" (set BUILD_DIR=%~2)
IF "%~1"=="-r" (set RELEASE_DIR=%~2)
IF "%~1"=="-g" (set NO_GUI="true")
SHIFT
GOTO parse
:endparse

call %SRC_DIR%\maracluster\admin\builders\setup_env.bat 32bit

set VCTARGET=%PROGRAM_FILES_DIR%\MSBuild\Microsoft.Cpp\v4.0\V%MSVC_VER%0

::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::: START INSTALL DEPENDENCIES ::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::

setlocal
set INSTALL_DIR=%BUILD_DIR%\tools
if not exist "%INSTALL_DIR%" (md "%INSTALL_DIR%")
if not exist "%RELEASE_DIR%" (md "%RELEASE_DIR%")

set ZIP_URL=http://downloads.sourceforge.net/sevenzip/7z920.exe
if not exist "%INSTALL_DIR%\7zip" (
  echo Downloading and installing 7-Zip
  call :downloadfile %ZIP_URL% %INSTALL_DIR%\7zip.exe
  "%INSTALL_DIR%\7zip.exe" /S /D=%INSTALL_DIR%\7zip
)
set ZIP_EXE="%INSTALL_DIR%\7zip\7z.exe"

set CMAKE_BASE=cmake-3.5.2-win32-x86
set CMAKE_URL=https://cmake.org/files/v3.5/%CMAKE_BASE%.zip
if not exist "%INSTALL_DIR%\%CMAKE_BASE%" (
  echo Downloading and installing CMake
  call :downloadfile %CMAKE_URL% %INSTALL_DIR%\cmake.zip
  %ZIP_EXE% x "%INSTALL_DIR%\cmake.zip" -o"%INSTALL_DIR%" -aoa -xr!doc > NUL
)
set CMAKE_EXE="%INSTALL_DIR%\%CMAKE_BASE%\bin\cmake.exe"

::: https://teamcity.labkey.org/viewType.html?buildTypeId=bt81 :::
::: without-t = without tests :::
set PWIZ_BASE=pwiz-src-without-t-3_0_19025_7f0e41d
set PWIZ_URL=https://teamcity.labkey.org/guestAuth/repository/download/bt81/691783:id/%PWIZ_BASE%.zip
set PWIZ_DIR=%INSTALL_DIR%\proteowizard
if not exist "%PWIZ_DIR%\lib" (
  echo Downloading and installing ProteoWizard
  if not exist "%PWIZ_DIR%" (
    call :downloadfile %PWIZ_URL% %INSTALL_DIR%\pwiz.tar.bz2
    if not exist "%INSTALL_DIR%\pwiz.tar" (
      %ZIP_EXE% x "%INSTALL_DIR%\pwiz.tar.bz2" -o"%INSTALL_DIR%" -aoa > NUL
    )
    %ZIP_EXE% x "%INSTALL_DIR%\pwiz.tar" -o"%PWIZ_DIR%" -aoa > NUL
  )
  cd /D "%PWIZ_DIR%"
  
  call quickbuild.bat -j4 --toolset=msvc-%MSVC_VER%.0 --i-agree-to-the-vendor-licenses ^
                pwiz/data/common//pwiz_data_common ^
                pwiz/data/identdata//pwiz_data_identdata ^
                pwiz/data/identdata//pwiz_data_identdata_version ^
                pwiz/data/msdata//pwiz_data_msdata ^
                pwiz/data/msdata//pwiz_data_msdata_version ^
                pwiz/data/msdata/mz5//pwiz_data_msdata_mz5 ^
                pwiz/data/proteome//pwiz_data_proteome ^
                pwiz/utility/chemistry//pwiz_utility_chemistry ^
                pwiz/utility/minimxml//pwiz_utility_minimxml ^
                pwiz/utility/misc//SHA1 ^
                pwiz/utility/misc//pwiz_utility_misc ^
                pwiz/data/vendor_readers//pwiz_data_vendor_readers ^
                pwiz/data/vendor_readers/Thermo//pwiz_reader_thermo ^
                pwiz_aux/msrc/utility/vendor_api/thermo//pwiz_vendor_api_thermo ^
                pwiz/data/vendor_readers/Shimadzu//pwiz_reader_shimadzu ^
                pwiz_aux/msrc/utility/vendor_api/Shimadzu//pwiz_vendor_api_shimadzu ^
                pwiz/data/vendor_readers/UIMF//pwiz_reader_uimf ^
                pwiz_aux/msrc/utility/vendor_api/UIMF//pwiz_vendor_api_uimf ^
                pwiz/data/vendor_readers/Agilent//pwiz_reader_agilent ^
                pwiz_aux/msrc/utility/vendor_api/Agilent//pwiz_vendor_api_agilent ^
                pwiz/data/vendor_readers/Waters//pwiz_reader_waters ^
                pwiz/data/vendor_readers/Bruker//pwiz_reader_bruker ^
                pwiz_aux/msrc/utility/vendor_api/Bruker//pwiz_vendor_api_bruker ^
                pwiz/data/vendor_readers/ABI/T2D//pwiz_reader_abi_t2d ^
                pwiz_aux/msrc/utility/vendor_api/ABI/T2D//pwiz_vendor_api_abi_t2d ^
                pwiz/data/vendor_readers/ABI//pwiz_reader_abi ^
                pwiz_aux/msrc/utility/vendor_api/ABI//pwiz_vendor_api_abi ^
                /ext/zlib//z ^
                /ext/hdf5//hdf5pp ^
                /ext/hdf5//hdf5 ^
                libraries/SQLite//sqlite3 ^
                libraries/SQLite//sqlite3pp ^
                /ext/boost//system ^
                /ext/boost//thread ^
                /ext/boost//chrono ^
                /ext/boost//regex ^
                /ext/boost//filesystem ^
                /ext/boost//iostreams ^
                /ext/boost//program_options ^
                /ext/boost//nowide ^
                /ext/boost//serialization > pwiz_installation.log 2>&1
  
  mkdir lib
  for /r build-nt-x86 %%x in (*.lib) do copy "%%x" lib\ /Y
  cd lib
  PowerShell "(Dir | Rename-Item -NewName { $_.Name -replace '-vc%MSVC_VER%0-mt','' })"
  PowerShell "(Dir | Rename-Item -NewName { $_.Name -replace 'libpwiz_','pwiz_' })"
  Ren libzlib.lib zlib.lib
  Ren libhdf5pp.lib hdf5pp.lib
  Ren libhdf5.lib hdf5.lib
  Ren libSHA1.lib SHA1.lib
  Ren libsqlite3pp.lib sqlite3pp.lib
  Ren libsqlite3.lib sqlite3.lib
  COPY ..\pwiz_aux\msrc\utility\vendor_api\Waters\vc12_x86\* .
  COPY ..\pwiz_aux\msrc\utility\vendor_api\Bruker\x86\baf2sql_c.* .
  cd ..

  mkdir include
  for /r pwiz %%x in (*.hpp, *.h) do copy "%%x" include\ /Y
)

set QT_BASE=qtbase-everywhere-src-5.11.2
set QT_URL=http://download.qt.io/official_releases/qt/5.11/5.11.2/submodules/%QT_BASE%.zip
set QT_DIR=%INSTALL_DIR%\%QT_BASE%
setlocal
set PATH=%PATH%;%QT_DIR%\bin
::: use multiple cores for nmake :::
set CL=/MP
if not "%NO_GUI%" == "true" (
  if not exist "%INSTALL_DIR%\Qt-dynamic" (
    echo Downloading Qt base
    call :downloadfile %QT_URL% %INSTALL_DIR%\qt.zip
    %ZIP_EXE% x "%INSTALL_DIR%\qt.zip" -o"%INSTALL_DIR%" -aoa
    
    cd /D "%QT_DIR%"

    configure -prefix "%INSTALL_DIR%\Qt-dynamic" -opensource -confirm-license -nomake tools -nomake examples -nomake tests
    
    echo Building Qt base, this may take some time..
    
    cd /D "%QT_DIR%"
    call :downloadfile http://code.qt.io/cgit/qt/qt5.git/plain/gnuwin32/bin/flex.exe %QT_DIR%\bin\flex.exe
    
    nmake
    nmake install
  )
)

::: Needed for CPack :::
set NSIS_DIR=%INSTALL_DIR%\nsis
set NSIS_URL=https://sourceforge.net/projects/nsis/files/NSIS%203/3.04/nsis-3.04-setup.exe/download
if not exist "%NSIS_DIR%" (
  echo Downloading and installing NSIS installer
  call :downloadfile %NSIS_URL% %INSTALL_DIR%\nsis.exe
  "%INSTALL_DIR%\nsis.exe" /S /D=%INSTALL_DIR%\nsis
)
set PATH=%PATH%;%INSTALL_DIR%\nsis

::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::: END INSTALL DEPENDENCIES ::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::



:::::::::::::::::::::::::::::::::::::::::
:::::::::::: START BUILD ::::::::::::::::
:::::::::::::::::::::::::::::::::::::::::

if not exist "%BUILD_DIR%" (md "%BUILD_DIR%")

::::::: Building maracluster :::::::
if not exist "%BUILD_DIR%\maracluster" (md "%BUILD_DIR%\maracluster")
cd /D "%BUILD_DIR%\maracluster....."
echo cmake maracluster.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -DBOOST_ROOT="%PWIZ_DIR%\libraries\boost_1_67_0" -DZLIB_INCLUDE_DIR="%PWIZ_DIR%\libraries\zlib-1.2.3" -DCMAKE_PREFIX_PATH="%PWIZ_DIR%" "%SRC_DIR%\maracluster"

echo build maracluster (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::::::: Building maracluster with vendor support :::::::
if not exist "%BUILD_DIR%\maracluster-vendor-support" (md "%BUILD_DIR%\maracluster-vendor-support")
cd /D "%BUILD_DIR%\maracluster-vendor-support"
echo cmake maracluster with vendor support.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -DBOOST_ROOT="%PWIZ_DIR%\libraries\boost_1_67_0" -DZLIB_INCLUDE_DIR="%PWIZ_DIR%\libraries\zlib-1.2.3" -DCMAKE_PREFIX_PATH="%PWIZ_DIR%" -DVENDOR_SUPPORT=ON "%SRC_DIR%\maracluster"

echo build maracluster with vendor support (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

if not "%NO_GUI%" == "true" (
  ::::::: Building maracluster with GUI :::::::
  if not exist "%BUILD_DIR%\maracluster-gui" (md "%BUILD_DIR%\maracluster-gui")
  cd /D "%BUILD_DIR%\maracluster-gui"
  echo cmake maracluster gui.....
  %CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -DBOOST_ROOT="%PWIZ_DIR%\libraries\boost_1_67_0" -DZLIB_INCLUDE_DIR="%PWIZ_DIR%\libraries\zlib-1.2.3" -DCMAKE_PREFIX_PATH="%PWIZ_DIR%;%INSTALL_DIR\Qt-dynamic" -DVENDOR_SUPPORT=OFF "%SRC_DIR%\maracluster\src\qt-gui"

  echo build maracluster gui (this will take a few minutes).....
  msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

  ::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
  ::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
)

:::::::::::::::::::::::::::::::::::::::
:::::::::::: END BUILD ::::::::::::::::
:::::::::::::::::::::::::::::::::::::::

echo Copying installers to %RELEASE_DIR%
copy "%BUILD_DIR%\maracluster\mar*.exe" "%RELEASE_DIR%"
copy "%BUILD_DIR%\maracluster-vendor-support\mar*.exe" "%RELEASE_DIR%"
copy "%BUILD_DIR%\maracluster-gui\mar*.exe" "%RELEASE_DIR%"

echo Finished buildscript execution in build directory %BUILD_DIR%

cd "%SRC_DIR%"

EXIT /B %errorlevel%


:downloadfile
PowerShell "[Net.ServicePointManager]::SecurityProtocol = 'tls12, tls11, tls'; (new-object System.Net.WebClient).DownloadFile('%1','%2')"
EXIT /B
