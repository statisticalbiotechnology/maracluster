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
set NO_GUI=false

:parse
IF "%~1"=="" GOTO endparse
IF "%~1"=="-s" (set SRC_DIR=%~2)
IF "%~1"=="-b" (set BUILD_DIR=%~2)
IF "%~1"=="-r" (set RELEASE_DIR=%~2)
IF "%~1"=="-g" (set NO_GUI=true)
SHIFT
GOTO parse
:endparse

del "%BUILD_DIR%\maracluster\mar*.exe" >nul 2>&1
del "%BUILD_DIR%\maracluster-gui\mar*.exe" >nul 2>&1

call %SRC_DIR%\maracluster\admin\builders\_init_msvc_.bat 32bit
if %ERRORLEVEL% NEQ 0 (
  EXIT /B %ERRORLEVEL%
)

::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::: START INSTALL DEPENDENCIES ::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::

setlocal

call %SRC_DIR%\maracluster\admin\builders\_urls_and_file_names_.bat

set INSTALL_DIR=%BUILD_DIR%\tools
if not exist "%INSTALL_DIR%" (md "%INSTALL_DIR%")
if not exist "%RELEASE_DIR%" (md "%RELEASE_DIR%")

set ZIP_DIR=%INSTALL_DIR%\%ZIP_BASE%
if not exist "%ZIP_DIR%" (
  echo Downloading and installing 7-Zip
  call :downloadfile %ZIP_URL% %INSTALL_DIR%\7zip.exe
  "%INSTALL_DIR%\7zip.exe" /S /D=%INSTALL_DIR%\7zip
)
set ZIP_EXE="%ZIP_DIR%\7z.exe"

set CMAKE_DIR=%INSTALL_DIR%\%CMAKE_BASE%
if not exist "%CMAKE_DIR%" (
  echo Downloading and installing CMake
  call :downloadfile %CMAKE_URL% %INSTALL_DIR%\cmake.zip
  %ZIP_EXE% x "%INSTALL_DIR%\cmake.zip" -o"%INSTALL_DIR%" -aoa -xr!doc > NUL
)
set CMAKE_EXE="%CMAKE_DIR%\bin\cmake.exe"

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
                pwiz/data/vendor_readers/UNIFI//pwiz_reader_unifi ^
                pwiz_aux/msrc/utility/vendor_api/UNIFI//pwiz_vendor_api_unifi ^
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
  for /r build-nt-x86 %%x in (*.lib) do copy "%%x" lib\ /Y > NUL
  cd lib
  PowerShell "(Dir | Rename-Item -NewName { $_.Name -replace '-vc%MSVC_VER%0-mt','' })"
  PowerShell "(Dir | Rename-Item -NewName { $_.Name -replace 'libpwiz_','pwiz_' })"
  Ren libzlib.lib zlib.lib
  Ren libhdf5pp.lib hdf5pp.lib
  Ren libhdf5.lib hdf5.lib
  Ren libSHA1.lib SHA1.lib
  Ren libsqlite3pp.lib sqlite3pp.lib
  Ren libsqlite3.lib sqlite3.lib
  COPY ..\pwiz_aux\msrc\utility\vendor_api\Waters\vc12_x86\* . > NUL
  COPY ..\pwiz_aux\msrc\utility\vendor_api\Bruker\x86\baf2sql_c.* . > NUL
  COPY ..\pwiz_aux\msrc\utility\vendor_api\Bruker\x86\timsdata.* . > NUL
  cd ..
  
  ::: Generate lib from dll for cdt.dll
  setlocal enableDelayedExpansion
  set DLL_BASE=cdt
  set DEF_FILE=!DLL_BASE!.def
  set write=0
  echo EXPORTS> "!DEF_FILE!"
  for /f "usebackq tokens=4" %%i in (`dumpbin /exports "!DLL_BASE!.dll"`) do if "!write!"=="1" (echo %%i >> "!DEF_FILE!") else (if %%i==name set write=1)
  lib /DEF:"!DEF_FILE!" /MACHINE:x86
  endlocal

  mkdir include
  for /r pwiz %%x in (*.hpp, *.h) do copy "%%x" include\ /Y > NUL
  
  ::: copy the boost::asio library, which is not included by the ProteoWizard boost tar but is needed for maracluster
  call :downloadfile "%BOOST_ASIO_URL%" %INSTALL_DIR%\boost_asio.zip
  %ZIP_EXE% x "%INSTALL_DIR%\boost_asio.zip" -o"%INSTALL_DIR%" -aoa > NUL
  PowerShell "Copy-Item -Path '%INSTALL_DIR%\%BOOST_ASIO_BASE%\boost' -Destination '%PWIZ_DIR%\libraries\boost_1_67_0' -Recurse -Force"
)

set QT_DIR=%INSTALL_DIR%\%QT_BASE%
if not "%NO_GUI%" == "true" (
  if not exist "%INSTALL_DIR%\Qt-dynamic" (
    ::: use multiple cores with jom instead of single-core nmake :::
    echo Downloading Jom
    call :downloadfile %JOM_URL% %INSTALL_DIR%\jom.zip
    %ZIP_EXE% x "%INSTALL_DIR%\jom.zip" -o"%INSTALL_DIR%\jom" -aoa > NUL
    
    echo Downloading Qt base
    call :downloadfile %QT_URL% %INSTALL_DIR%\qt.zip
    %ZIP_EXE% x "%INSTALL_DIR%\qt.zip" -o"%INSTALL_DIR%" -aoa > NUL
    
    cd /D "%QT_DIR%"

    call configure.bat -prefix "%INSTALL_DIR%\Qt-dynamic" -opensource -confirm-license -nomake tools -nomake examples -nomake tests -release > qt_config.log 2>&1
    
    echo Building Qt base, this may take some time..
    
    cd /D "%QT_DIR%"
    %ZIP_EXE% x "%SRC_DIR%\maracluster\admin\gnuwin32_bin.zip" -o"%INSTALL_DIR%" -aoa > NUL
    setlocal
    set "PATH=%PATH%;%INSTALL_DIR%\jom;%QT_DIR%\bin;%INSTALL_DIR%\gnuwin32_bin"
    jom > qt_installation.log 2>&1
    jom install
    endlocal
  )
)

::: Needed for CPack :::
set NSIS_DIR=%INSTALL_DIR%\nsis
if not exist "%NSIS_DIR%" (
  echo Downloading and installing NSIS installer
  call :downloadfile "%NSIS_URL%" %INSTALL_DIR%\nsis.exe
  "%INSTALL_DIR%\nsis.exe" /S /D=%INSTALL_DIR%\nsis
)
setlocal
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

echo build maracluster.....
msbuild PACKAGE.vcxproj /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:Configuration=%BUILD_TYPE% /m

if not "%NO_GUI%" == "true" (
  ::::::: Building maracluster with GUI :::::::
  if not exist "%BUILD_DIR%\maracluster-gui" (md "%BUILD_DIR%\maracluster-gui")
  cd /D "%BUILD_DIR%\maracluster-gui"
  echo cmake maracluster gui.....
  %CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -DBOOST_ROOT="%PWIZ_DIR%\libraries\boost_1_67_0" -DZLIB_INCLUDE_DIR="%PWIZ_DIR%\libraries\zlib-1.2.3" -DCMAKE_PREFIX_PATH="%PWIZ_DIR%;%INSTALL_DIR%\Qt-dynamic" -DVENDOR_SUPPORT=OFF "%SRC_DIR%\maracluster\src\qt-gui"

  echo build maracluster gui.....
  msbuild PACKAGE.vcxproj /p:Configuration=%BUILD_TYPE% /m

  ::msbuild INSTALL.vcxproj /p:Configuration=%BUILD_TYPE% /m
  ::msbuild RUN_TESTS.vcxproj /p:Configuration=%BUILD_TYPE% /m
)

echo Finished buildscript execution in build directory %BUILD_DIR%

:::::::::::::::::::::::::::::::::::::::
:::::::::::: END BUILD ::::::::::::::::
:::::::::::::::::::::::::::::::::::::::

echo Copying installers to %RELEASE_DIR%
set /A exit_code=0
call :copytorelease "%BUILD_DIR%\maracluster\mar*.exe"

if not "%NO_GUI%" == "true" (
  call :copytorelease "%BUILD_DIR%\maracluster-gui\mar*.exe"
)

cd /D "%SRC_DIR%"

EXIT /B %exit_code%

:downloadfile
echo Downloading %1 to %2
PowerShell "[Net.ServicePointManager]::SecurityProtocol = 'tls12, tls11, tls'; (new-object System.Net.WebClient).DownloadFile('%1','%2')"
EXIT /B

:copytorelease
echo Copying "%1" to "%RELEASE_DIR%"
xcopy %1 "%RELEASE_DIR%" /Y
dir %1 /b /a-d >nul 2>&1
set /A exit_code=exit_code+%ERRORLEVEL%
EXIT /B
