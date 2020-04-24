::: Centralized place for urls and files for all windows builders ...
::: please do not change compression type in urls, since decompression is
::: hardcoded in the respective buiding scripts

::: 7-zip
set ZIP_BASE=7zi
set ZIP_URL=https://www.7-zip.org/a/7z1900.exe

::: CMake
set CMAKE_VERSION=3.16.6
set CMAKE_BASE=cmake-%CMAKE_VERSION%-win32-x86
set CMAKE_URL=https://github.com/Kitware/CMake/releases/download/v%CMAKE_VERSION%/%CMAKE_BASE%.zip

::: Proteowizard
::: https://teamcity.labkey.org/viewType.html?buildTypeId=bt81 :::
::: without-t = without tests :::
set PWIZ_VERSION_URL=https://teamcity.labkey.org/guestAuth/repository/download/bt81/.lastSuccessful/VERSION
call :downloadfile %PWIZ_VERSION_URL% %INSTALL_DIR%\VERSION
set /p PWIZ_VERSION_STRING=<%INSTALL_DIR%\VERSION
set PWIZ_BASE=pwiz-src-without-t-%PWIZ_VERSION_STRING: =_%
set PWIZ_URL=https://teamcity.labkey.org/guestAuth/repository/download/bt81/.lastSuccessful/%PWIZ_BASE%.tar.bz2

::: Qt
set QT_BASE=qtbase-opensource-src-5.9.9
set QT_URL=https://download.qt.io/archive/qt/5.9/5.9.9/submodules/%QT_BASE%.zip
set JOM_URL=http://download.qt.io/official_releases/jom/jom_1_1_3.zip

::: NSIS
set NSIS_URL=https://sourceforge.net/projects/nsis/files/NSIS 3/3.04/nsis-3.04-setup.exe/download

EXIT /B

:downloadfile
PowerShell "[Net.ServicePointManager]::SecurityProtocol = 'tls12, tls11, tls'; (new-object System.Net.WebClient).DownloadFile('%1','%2')"
EXIT /B

