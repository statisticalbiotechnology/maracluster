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

call download_file_macro.bat >NUL

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
  "%INSTALL_DIR%\7zip.exe" /S /D=%ZIP_DIR%
)
set ZIP_EXE="%ZIP_DIR%\7z.exe"

set CMAKE_DIR=%INSTALL_DIR%\%CMAKE_BASE%
if not exist "%CMAKE_DIR%" (
  echo Downloading and installing CMake
  call :downloadfile %CMAKE_URL% %INSTALL_DIR%\cmake.zip
  %ZIP_EXE% x "%INSTALL_DIR%\cmake.zip" -o"%INSTALL_DIR%" -aoa -xr!doc > NUL
)
set CMAKE_EXE="%CMAKE_DIR%\bin\cmake.exe"

call %SRC_DIR%\maracluster\admin\builders\install_proteowizard.bat 32bit

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
%CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -DBOOST_ROOT="%PWIZ_DIR%\libraries\%BOOST_BASE%" -DZLIB_INCLUDE_DIR="%PWIZ_DIR%\libraries\%ZLIB_BASE%" -DCMAKE_PREFIX_PATH="%PWIZ_DIR%" "%SRC_DIR%\maracluster"

echo build maracluster.....
msbuild PACKAGE.vcxproj /p:Configuration=%BUILD_TYPE% /m

::msbuild INSTALL.vcxproj /p:Configuration=%BUILD_TYPE% /m
::msbuild RUN_TESTS.vcxproj /p:Configuration=%BUILD_TYPE% /m

if not "%NO_GUI%" == "true" (
  ::::::: Building maracluster with GUI :::::::
  if not exist "%BUILD_DIR%\maracluster-gui" (md "%BUILD_DIR%\maracluster-gui")
  cd /D "%BUILD_DIR%\maracluster-gui"
  echo cmake maracluster gui.....
  %CMAKE_EXE% -G "Visual Studio %MSVC_VER%" -DBOOST_ROOT="%PWIZ_DIR%\libraries\%BOOST_BASE%" -DZLIB_INCLUDE_DIR="%PWIZ_DIR%\libraries\%ZLIB_BASE%" -DCMAKE_PREFIX_PATH="%PWIZ_DIR%;%INSTALL_DIR%\Qt-dynamic" -DVENDOR_SUPPORT=OFF "%SRC_DIR%\maracluster\src\qt-gui"

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


:copytorelease
echo Copying "%1" to "%RELEASE_DIR%"
xcopy %1 "%RELEASE_DIR%" /Y
dir %1 /b /a-d >nul 2>&1
set /A exit_code=exit_code+%ERRORLEVEL%
EXIT /B
