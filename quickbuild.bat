@echo off
set MSVC_VER=12
set VCTARGET=C:\Program Files\MSBuild\Microsoft.Cpp\v4.0\V%MSVC_VER%0
set SRC_DIR=%~dp0
set BUILD_DIR=%SRC_DIR%\..\build
set RELEASE_DIR=%SRC_DIR%\..\release
set BUILD_TYPE=Release

:parse
IF "%~1"=="" GOTO endparse
IF "%~1"=="-s" (set SRC_DIR=%~2)
IF "%~1"=="-b" (set BUILD_DIR=%~2)
IF "%~1"=="-r" (set RELEASE_DIR=%~2)
SHIFT
GOTO parse
:endparse

:: use the VS command prompt settings to set-up paths for compiler and builder
call "C:\Program Files\Microsoft Visual Studio %MSVC_VER%.0\Common7\Tools\VsDevCmd.bat"

::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::::::::::: START INSTALL DEPENDENCIES ::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::::::::::

setlocal
set INSTALL_DIR=%BUILD_DIR%\tools
if not exist "%INSTALL_DIR%" (md "%INSTALL_DIR%")
if not exist "%RELEASE_DIR%" (md "%RELEASE_DIR%")

set CMAKE_URL=http://www.cmake.org/files/v2.8/cmake-2.8.12.1-win32-x86.exe
if not exist "%INSTALL_DIR%\cmake" (
  echo Downloading and installing CMake
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%CMAKE_URL%','%INSTALL_DIR%\cmake.exe')"
  "%INSTALL_DIR%\cmake.exe" /S /D=%INSTALL_DIR%\cmake
)
set CMAKE_EXE="%INSTALL_DIR%\cmake\bin\cmake.exe"

set ZIP_URL=http://downloads.sourceforge.net/sevenzip/7z920.exe
if not exist "%INSTALL_DIR%\7zip" (
  echo Downloading and installing 7-Zip
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%ZIP_URL%','%INSTALL_DIR%\7zip.exe')"
  "%INSTALL_DIR%\7zip.exe" /S /D=%INSTALL_DIR%\7zip
)
set ZIP_EXE="%INSTALL_DIR%\7zip\7z.exe"

::: Needed for proteowizard
set SVN_URL="https://sourceforge.net/projects/win32svn/files/latest/download"
WHERE svn >nul 2>&1
IF ERRORLEVEL 1 (
  echo Downloading and installing SVN
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%SVN_URL%','%INSTALL_DIR%\svn.msi')"
  cd /D "%INSTALL_DIR%"
  msiexec /i svn.msi /quiet /Li svn_install.log
)

set PWIZ_DIR=%INSTALL_DIR%\proteowizard
if not exist "%PWIZ_DIR%" (
  echo Downloading and installing ProteoWizard
  REV="7692"
  svn co -r %REV% --depth immediates https://svn.code.sf.net/p/proteowizard/code/trunk/pwiz %PWIZ_DIR%
  svn update -r %REV% --set-depth infinity %PWIZ_DIR%\pwiz
  svn update -r %REV% --set-depth infinity %PWIZ_DIR%\pwiz_aux
  svn update -r %REV% --set-depth infinity %PWIZ_DIR%\libraries
  cd /D %PWIZ_DIR%
  call quickbuild.bat --i-agree-to-the-vendor-licenses ^
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
                pwiz/data/vendor_readers/UIMF//pwiz_reader_uimf ^
                pwiz/data/vendor_readers/Agilent//pwiz_reader_agilent ^
                pwiz/data/vendor_readers/Waters//pwiz_reader_waters ^
                pwiz/data/vendor_readers/Bruker//pwiz_reader_bruker ^
                pwiz/data/vendor_readers/ABI/T2D//pwiz_reader_abi_t2d ^
                pwiz/data/vendor_readers/ABI//pwiz_reader_abi ^
                /ext/zlib//z ^
                /ext/hdf5//hdf5pp ^
                /ext/hdf5//hdf5 ^
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
  PowerShell "(Dir | Rename-Item -NewName { $_.Name -replace '-vc120-mt','' })"
  PowerShell "(Dir | Rename-Item -NewName { $_.Name -replace 'libpwiz_','pwiz_' })"
  Ren libzlib.lib zlib.lib
  Ren libhdf5pp.lib hdf5pp.lib
  Ren libhdf5.lib hdf5.lib
  Ren libSHA1.lib SHA1.lib
  cd ..

  mkdir include
  for /r pwiz %%x in (*.hpp, *.h) do copy "%%x" include\ /Y
)

::: Needed for CPack :::
set NSIS_DIR=%INSTALL_DIR%\nsis
set NSIS_URL=http://downloads.sourceforge.net/project/nsis/NSIS 3 Pre-release/3.0a2/nsis-3.0a2-setup.exe
if not exist "%NSIS_DIR%" (
  echo Downloading and installing NSIS installer
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%NSIS_URL%','%INSTALL_DIR%\nsis.exe')"
  "%INSTALL_DIR%\nsis.exe" /S /D=%INSTALL_DIR%\nsis
)

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
%CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -DBOOST_ROOT="%PWIZ_DIR%\libraries\boost_1_56_0" -DZLIB_INCLUDE_DIR="%PWIZ_DIR%\libraries\zlib-1.2.3" -DCMAKE_PREFIX_PATH="%PWIZ_DIR%" "%SRC_DIR%"

echo build maracluster (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

:::::::::::::::::::::::::::::::::::::::
:::::::::::: END BUILD ::::::::::::::::
:::::::::::::::::::::::::::::::::::::::

echo Copying installers to %RELEASE_DIR%
copy "%BUILD_DIR%\maracluster\mar*.exe" "%RELEASE_DIR%"

echo Finished buildscript execution in build directory %BUILD_DIR%

cd "%SRC_DIR%"
