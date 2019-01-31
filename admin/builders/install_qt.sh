tools_dir=$1

# change directory to the build director
cd ${tools_dir}

linux_qt=qtbase-everywhere-src-5.11.2
if [ ! -d Qt-dynamic ]; then

  if [ ! -f ${linux_qt}.tar.xz ]; then
    wget http://download.qt.io/official_releases/qt/5.11/5.11.2/submodules/${linux_qt}.tar.xz
  fi

  tar xf ${linux_qt}.tar.xz
  cd ${linux_qt}

  ./configure -prefix ${tools_dir}/Qt-dynamic -opensource -confirm-license -nomake tools -nomake examples -nomake tests
  
  echo "Building Qt base, this may take some time.."
  
  make -j4 > ../qt_installation.log 2>&1
  make install -j4
fi

cd ${tools_dir}

# Qt5 requires CMake >= 3.5
function version_lt() { test "$(echo "$@" | tr " " "\n" | sort -rV | head -n 1)" != "$1"; }

if version_lt $(cmake --version | head -n1 | cut -f3 -d ' ') "3.5"; then
  wget https://github.com/Kitware/CMake/releases/download/v3.13.3/cmake-3.13.3-Linux-x86_64.sh
  bash cmake-3.13.3-Linux-x86_64.sh --skip-license --exclude-subdir
  CMAKE_BINARY=${build_dir}/tools/bin/cmake
else
  CMAKE_BINARY=cmake 
fi

export CMAKE_BINARY