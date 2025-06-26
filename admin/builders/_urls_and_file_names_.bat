::: Centralized place for urls and files for all windows builders ...
::: please do not change compression type in urls, since decompression is
::: hardcoded in the respective buiding scripts

::: 7-zip
set ZIP_BASE=7zip
set ZIP_URL=https://sourceforge.net/projects/sevenzip/files/7-Zip/24.09/7z2409.exe/download

::: CMake
set CMAKE_VERSION=3.23.1
set CMAKE_BASE=cmake-%CMAKE_VERSION%-windows-i386
set CMAKE_URL=https://github.com/Kitware/CMake/releases/download/v%CMAKE_VERSION%/%CMAKE_BASE%.zip

::: Proteowizard
::: https://teamcity.labkey.org/viewType.html?buildTypeId=bt81 :::
::: without-t = without tests :::
set PWIZ_VERSION_URL=https://proteowizard.sourceforge.io/releases/bt81.xml
call %~dp0\download_file.bat %PWIZ_VERSION_URL% bt81.xml
call :extractbuildinfo
echo ✅ PWIZ_BUILD_ID: %PWIZ_BUILD_ID%
echo ✅ PWIZ_FILE_NAME: %PWIZ_FILE_NAME%
set PWIZ_URL=https://mc-tca-01.s3.us-west-2.amazonaws.com/ProteoWizard/bt81/%PWIZ_BUILD_ID%/%PWIZ_FILE_NAME%

::: remember to update the boost and zlib version when ProteoWizard updates their versions
set BOOST_BASE=boost_1_86_0
set ZLIB_BASE=zlib-1.2.3

::: Boost asio library
set BOOST_ASIO_BASE=boost_asio_1_30_2
set BOOST_ASIO_URL=https://sourceforge.net/projects/asio/files/asio/1.30.2%20(Stable)/%BOOST_ASIO_BASE%.zip/download

::: Boost unordered library
set BOOST_UNORDERED_BASE=unordered-boost-1.86.0
set BOOST_UNORDERED_URL=https://github.com/boostorg/unordered/archive/refs/tags/boost-1.86.0.zip

::: Qt
set QT_BASE=qtbase-opensource-src-5.9.9
set QT_URL=https://download.qt.io/archive/qt/5.9/5.9.9/submodules/%QT_BASE%.zip
set JOM_URL=http://download.qt.io/official_releases/jom/jom_1_1_3.zip

::: NSIS
set NSIS_URL=https://sourceforge.net/projects/nsis/files/NSIS%203/3.11/nsis-3.11-setup.exe/download

EXIT /B


:: Macro to extract BUILD_ID and FILE_NAME using PowerShell regex
:extractbuildinfo
::: Step 1: Run PowerShell and output extracted strings to temp files
PowerShell "$content = Get-Content -Raw 'bt81.xml'; if ($content -match 'id:(\d+)/artifacts/content/(pwiz-src-without-t-[^<]+?\.tar\.bz2)') { $matches[1] | Out-File -Encoding ASCII 'build_id.tmp' }"
PowerShell "$content = Get-Content -Raw 'bt81.xml'; if ($content -match 'id:(\d+)/artifacts/content/(pwiz-src-without-t-[^<]+?\.tar\.bz2)') { $matches[2] | Out-File -Encoding ASCII 'file_name.tmp' }"
::: Step 2: Read from the temp files into batch variables
set /p PWIZ_BUILD_ID=<build_id.tmp
set /p PWIZ_FILE_NAME=<file_name.tmp
EXIT /B
