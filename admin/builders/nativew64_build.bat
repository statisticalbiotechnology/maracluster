:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
:::                                                                           :::
::: tested on a Windows 10 64-bit installation with VS2015 Community          :::
:::                                                                           :::
::: N.B.: VS2015 Express cannot compile x64 natively and will break the build :::
:::                                                                           :::
::: To support reading Thermo RAW Files you need the MSFileReader executable, :::
::: which can be downloaded from https://thermo.flexnetoperations.com/.       :::
::: Search for "MSFileReader" and download and install the 3.0 SP3 version    :::
:::                                                                           :::
::: To support reading Waters RAW Files, Visual C++ runtime is necessary to   :::
::: run http://www.microsoft.com/en-us/download/details.aspx?id=30679         :::
:::                                                                           :::
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

@echo off
set MSVC_VER=14
set VCTARGET=C:\Program Files (x86)\MSBuild\Microsoft.cpp\v4.0\V%MSVC_VER%0
set SRC_DIR=%~dp0..\..\..\
set BUILD_DIR=%SRC_DIR%\build\win64
set RELEASE_DIR=%SRC_DIR%\release\win64
set BUILD_TYPE=Release

:parse
IF "%~1"=="" GOTO endparse
IF "%~1"=="-s" (set SRC_DIR=%~2)
IF "%~1"=="-b" (set BUILD_DIR=%~2)
IF "%~1"=="-r" (set RELEASE_DIR=%~2)
SHIFT
GOTO parse
:endparse

:: del "%BUILD_DIR%\maracluster\mar*.exe"
:: del "%BUILD_DIR%\maracluster-vendor-support\mar*.exe"
:: del "%BUILD_DIR%\maracluster-gui\mar*.exe"

:: use the VS command prompt settings to set-up paths for compiler and builder
call "C:\Program Files (x86)\Microsoft Visual Studio %MSVC_VER%.0\Common7\Tools\VsDevCmd.bat"
call "C:\Program Files (x86)\Microsoft Visual Studio %MSVC_VER%.0\VC\vcvarsall.bat" amd64

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
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%ZIP_URL%','%INSTALL_DIR%\7zip.exe')"
  "%INSTALL_DIR%\7zip.exe" /S /D=%INSTALL_DIR%\7zip
)
set ZIP_EXE="%INSTALL_DIR%\7zip\7z.exe"

set CMAKE_BASE=cmake-3.5.2-win32-x86
set CMAKE_URL=https://cmake.org/files/v3.5/%CMAKE_BASE%.zip
if not exist "%INSTALL_DIR%\%CMAKE_BASE%" (
  echo Downloading and installing CMake
  PowerShell "(new-object System.Net.WebClient).DownloadFile('%CMAKE_URL%','%INSTALL_DIR%\cmake.zip')"
  %ZIP_EXE% x "%INSTALL_DIR%\cmake.zip" -o"%INSTALL_DIR%" -aoa -xr!doc > NUL
)
set CMAKE_EXE="%INSTALL_DIR%\%CMAKE_BASE%\bin\cmake.exe"

::: have not figured out how to install SVN in a custom folder, msiexec + INSTALLDIR/TARGETDIR doesn't work
set SVN_URL=https://sourceforge.net/projects/win32svn/files/latest/download
WHERE svn >nul 2>&1
IF ERRORLEVEL 1 (
  echo Downloading and installing SVN
  PowerShell "[Net.ServicePointManager]::SecurityProtocol = 'tls12, tls11, tls'; (new-object System.Net.WebClient).DownloadFile('%SVN_URL%','%INSTALL_DIR%\svn.msi')"
  cd /D "%INSTALL_DIR%"
  msiexec /i svn.msi /quiet /Li svn_install.log
)
setlocal
set PATH=%PATH%;C:\Program Files (x86)\Subversion\bin

set PWIZ_DIR=%INSTALL_DIR%\proteowizard
set REV=10210
set SHIMADZU_JAMFILE=%PWIZ_DIR%\pwiz\data\vendor_readers\Shimadzu\Jamfile.jam
set SHIMADZU_JAMFILE2=%PWIZ_DIR%\pwiz_aux\msrc\utility\vendor_api\Shimadzu\Jamfile.jam
set HDF5_FILE=%PWIZ_DIR%\libraries\hdf5-1.8.7\src\windows\H5pubconf.h
if not exist "%PWIZ_DIR%\lib" (
  echo Downloading and installing ProteoWizard
  if not exist "%PWIZ_DIR%" (
    svn co -r %REV% --depth immediates svn://svn.code.sf.net/p/proteowizard/code/trunk/pwiz %PWIZ_DIR%
    svn update -r %REV% --set-depth infinity %PWIZ_DIR%\pwiz
    svn update -r %REV% --set-depth infinity %PWIZ_DIR%\pwiz_aux
    svn update -r %REV% --set-depth infinity %PWIZ_DIR%\libraries
  )
  
  ::: There is an issue with some TIMEZONE variables in VS2015 with hdf5, solved in v1.8.16. :::
  if not exist "%PWIZ_DIR%\libraries\hdf5-1.8.7.tar" (
    %ZIP_EXE% x "%PWIZ_DIR%\libraries\hdf5-1.8.7.tar.bz2" -o"%PWIZ_DIR%\libraries"
  )
  if not exist "%PWIZ_DIR%\libraries\hdf5-1.8.7" (
    %ZIP_EXE% x "%PWIZ_DIR%\libraries\hdf5-1.8.7.tar" -o"%PWIZ_DIR%\libraries"
    PowerShell "(Get-Content '%HDF5_FILE%') | ForEach-Object { $_ -replace '#define H5_HAVE_TIMEZONE 1', '/* #define H5_HAVE_TIMEZONE 1 */' } | Set-Content '%HDF5_FILE%'"
  )
  
  ::: The Shimadzu API build has to be changed to include VS2015 in the Jamfile :::
  PowerShell "(Get-Content '%SHIMADZU_JAMFILE%') | ForEach-Object { $_ -replace '12.0\" \"12.0express', '12.0\" \"14.0\" \"14.0express\" \"12.0express' } | Set-Content '%SHIMADZU_JAMFILE%'"
  PowerShell "(Get-Content '%SHIMADZU_JAMFILE2%') | ForEach-Object { $_ -replace '12.0\" \"12.0express', '12.0\" \"14.0\" \"14.0express\" \"12.0express' } | Set-Content '%SHIMADZU_JAMFILE2%'"
  
  cd /D "%PWIZ_DIR%"
  echo Starting Proteowizard and Boost build, this can take a while...
  call quickbuild.bat address-model=64 -j4 --toolset=msvc-%MSVC_VER%.0 --i-agree-to-the-vendor-licenses ^
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
                /ext/hdf5//hdf5 ^
                /ext/hdf5//hdf5pp ^
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
  for /r build-nt-x86_64 %%x in (*.lib) do copy "%%x" lib\ /Y
  cd lib
  PowerShell "(Dir | Rename-Item -NewName { $_.Name -replace '-vc%MSVC_VER%0-mt','' })"
  PowerShell "(Dir | Rename-Item -NewName { $_.Name -replace 'libpwiz_','pwiz_' })"
  Ren libzlib.lib zlib.lib
  Ren libhdf5pp.lib hdf5pp.lib
  Ren libhdf5.lib hdf5.lib
  Ren libSHA1.lib SHA1.lib
  Ren libsqlite3pp.lib sqlite3pp.lib
  Ren libsqlite3.lib sqlite3.lib
  ::: these DLLs might not work, as they are for VS2013 :::
  COPY ..\pwiz_aux\msrc\utility\vendor_api\Waters\vc12_x64\* .
  COPY ..\pwiz_aux\msrc\utility\vendor_api\Bruker\x64\baf2sql_c.* .
  cd ..

  mkdir include
  for /r pwiz %%x in (*.hpp, *.h) do copy "%%x" include\ /Y
)

::: Needed for CPack :::
set NSIS_DIR=%INSTALL_DIR%\nsis
set NSIS_URL=https://sourceforge.net/projects/nsis/files/NSIS 3 Pre-release/3.0rc1/nsis-3.0rc1-setup.exe/download
if not exist "%NSIS_DIR%" (
  echo Downloading and installing NSIS installer
  PowerShell "[Net.ServicePointManager]::SecurityProtocol = 'tls12, tls11, tls'; (new-object System.Net.WebClient).DownloadFile('%NSIS_URL%','%INSTALL_DIR%\nsis.exe')"
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
cd /D "%BUILD_DIR%\maracluster"
echo cmake maracluster.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER% Win64" -DBOOST_ROOT="%PWIZ_DIR%\libraries\boost_1_56_0" -DZLIB_INCLUDE_DIR="%PWIZ_DIR%\libraries\zlib-1.2.3" -DCMAKE_PREFIX_PATH="%PWIZ_DIR%" "%SRC_DIR%\maracluster"

echo build maracluster (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::::::: Building maracluster with vendor support :::::::
if not exist "%BUILD_DIR%\maracluster-vendor-support" (md "%BUILD_DIR%\maracluster-vendor-support")
cd /D "%BUILD_DIR%\maracluster-vendor-support"
echo cmake maracluster with vendor support.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER% Win64" -DBOOST_ROOT="%PWIZ_DIR%\libraries\boost_1_56_0" -DZLIB_INCLUDE_DIR="%PWIZ_DIR%\libraries\zlib-1.2.3" -DCMAKE_PREFIX_PATH="%PWIZ_DIR%" -DVENDOR_SUPPORT=ON "%SRC_DIR%\maracluster"

echo build maracluster with vendor support (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::::::: Building maracluster with GUI :::::::
if not exist "%BUILD_DIR%\maracluster-gui" (md "%BUILD_DIR%\maracluster-gui")
cd /D "%BUILD_DIR%\maracluster-gui"
echo cmake maracluster gui.....
%CMAKE_EXE% -G "Visual Studio %MSVC_VER% Win64" -DBOOST_ROOT="%PWIZ_DIR%\libraries\boost_1_56_0" -DZLIB_INCLUDE_DIR="%PWIZ_DIR%\libraries\zlib-1.2.3" -DCMAKE_PREFIX_PATH="%PWIZ_DIR%;C:\Qt\Qt5.11.2\5.11.2\msvc2015_64" -DVENDOR_SUPPORT=OFF "%SRC_DIR%\maracluster\src\qt-gui"

echo build maracluster gui (this will take a few minutes).....
msbuild PACKAGE.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:VCTargetsPath="%VCTARGET%" /p:Configuration=%BUILD_TYPE% /m

:::::::::::::::::::::::::::::::::::::::
:::::::::::: END BUILD ::::::::::::::::
:::::::::::::::::::::::::::::::::::::::

echo Copying installers to %RELEASE_DIR%
copy "%BUILD_DIR%\maracluster\mar*.exe" "%RELEASE_DIR%"
copy "%BUILD_DIR%\maracluster-vendor-support\mar*.exe" "%RELEASE_DIR%"
copy "%BUILD_DIR%\maracluster-gui\mar*.exe" "%RELEASE_DIR%"

echo Finished buildscript execution in build directory %BUILD_DIR%

cd "%SRC_DIR%"
