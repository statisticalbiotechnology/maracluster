set BUILD_TARGET=%1
if "%BUILD_TARGET%" == "32bit" (
  set ADDRESS_MODEL=32
  set TARGET_ARCH=x86
) else (
  set ADDRESS_MODEL=64
  set TARGET_ARCH=x64
)

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
  echo Starting Proteowizard and Boost build, this can take a while...
  call quickbuild.bat address-model=%ADDRESS_MODEL% -j4 --toolset=msvc-%MSVC_VER%.0 --i-agree-to-the-vendor-licenses ^
                pwiz/data/common//pwiz_data_common ^
                pwiz/data/identdata//pwiz_data_identdata ^
                pwiz/data/identdata//pwiz_data_identdata_version ^
                pwiz/data/msdata//pwiz_data_msdata_core ^
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
    
  echo Copying ProteoWizard libraries to lib folder
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
  
  ::: these DLLs might not work, as they are for VS2013 :::
  echo Copying vendor libraries to lib folder
  COPY ..\pwiz_aux\msrc\utility\vendor_api\Waters\vc12_%TARGET_ARCH%\* . > NUL
  COPY ..\pwiz_aux\msrc\utility\vendor_api\Bruker\%TARGET_ARCH%\baf2sql_c.* . > NUL
  COPY ..\pwiz_aux\msrc\utility\vendor_api\Bruker\%TARGET_ARCH%\timsdata.* . > NUL
  
  ::: Generate lib from dll for cdt.dll
  echo Generate lib from cdt.dll
  setlocal enableDelayedExpansion
  set DLL_BASE=cdt
  set DEF_FILE=!DLL_BASE!.def
  set write=0
  echo EXPORTS> "!DEF_FILE!"
  for /f "usebackq tokens=4" %%i in (`dumpbin /exports "!DLL_BASE!.dll"`) do if "!write!"=="1" (echo %%i >> "!DEF_FILE!") else (if %%i==name set write=1)
  lib /DEF:"!DEF_FILE!" /MACHINE:%TARGET_ARCH%
  endlocal
  
  cd ..
  
  echo Copying ProteoWizard header files to include folder
  mkdir include
  for /r pwiz %%x in (*.hpp, *.h) do copy "%%x" include\ /Y > NUL
  
  ::: copy the boost::asio library, which is not included by the ProteoWizard boost tar but is needed for maracluster
  call :downloadfile "%BOOST_ASIO_URL%" %INSTALL_DIR%\boost_asio.zip
  %ZIP_EXE% x "%INSTALL_DIR%\boost_asio.zip" -o"%INSTALL_DIR%" -aoa > NUL
  PowerShell "Copy-Item -Path '%INSTALL_DIR%\%BOOST_ASIO_BASE%\boost' -Destination '%PWIZ_DIR%\libraries\boost_1_86_0' -Recurse -Force"
)

EXIT /B


:downloadfile
PowerShell "[Net.ServicePointManager]::SecurityProtocol = 'tls12, tls11, tls'; (new-object System.Net.WebClient).DownloadFile('%1','%2')"
EXIT /B
