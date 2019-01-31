set BUILD_TARGET="%~1"
set MSVC_VER=0

:: use VS2015 if available
REG QUERY HKEY_CLASSES_ROOT\VisualStudio.DTE.14.0 > nul 2> nul
if %ERRORLEVEL% EQU 0 (
  echo Using Visual Studio 2015
  set MSVC_VER=14
) else (
  :: reset ERRORLEVEL to 0. N.B. set ERRORLEVEL=0 will permanently set it to 0
  cd .
)

:: fall back to VS2013 is available
if %MSVC_VER% EQU 0 (
  REG QUERY HKEY_CLASSES_ROOT\VisualStudio.DTE.12.0 > nul 2> nul
  if %ERRORLEVEL% EQU 0 (
    echo Using Visual Studio 2013
    set MSVC_VER=12
  ) else (
    :: reset ERRORLEVEL to 0
    cd .
  )
)

if %MSVC_VER% EQU 0 (
  echo Could not find a suitable Visual Studio version; supported versions: VS2013, VS2015
  EXIT /B 1  
)

set PROGRAM_FILES_DIR=C:\Program Files
set BUILD_PLATFORM=32bit
REG QUERY HKEY_LOCAL_MACHINE\SOFTWARE\WOW6432Node\Microsoft\VisualStudio\%MSVC_VER%.0\Setup\VS > nul 2> nul
if %ERRORLEVEL% EQU 0 (
  echo platform detected: 64-bit
  set BUILD_PLATFORM=64bit
  set "PROGRAM_FILES_DIR=C:\Program Files (x86)"
) else (
  :: reset ERRORLEVEL to 0
  cd .
)

:: use the VS command prompt settings to set-up paths for compiler and builder
:: see https://msdn.microsoft.com/en-us/library/f2ccy3wt.aspx for possible vcvarsall.bat arguments
if not defined DevEnvDir (
  call "%PROGRAM_FILES_DIR%\Microsoft Visual Studio %MSVC_VER%.0\Common7\Tools\VsDevCmd.bat"
  if "%BUILD_TARGET%" == "64bit" (
    if "%BUILD_PLATFORM%" == "64bit" (
      echo Setting variables for 64-bit
      call "%PROGRAM_FILES_DIR%\Microsoft Visual Studio %MSVC_VER%.0\VC\vcvarsall.bat" amd64
    ) else (
      echo Setting variables for 32-bit
      call "%PROGRAM_FILES_DIR%\Microsoft Visual Studio %MSVC_VER%.0\VC\vcvarsall.bat" x86_amd64
    )
  ) else (
    if "%BUILD_PLATFORM%" == "64bit" (
      echo Setting variables for 64-bit
      call "%PROGRAM_FILES_DIR%\Microsoft Visual Studio %MSVC_VER%.0\VC\vcvarsall.bat" amd64_x86
    ) else (
      echo Setting variables for 32-bit
      call "%PROGRAM_FILES_DIR%\Microsoft Visual Studio %MSVC_VER%.0\VC\vcvarsall.bat" x86
    )
  ) 
)
