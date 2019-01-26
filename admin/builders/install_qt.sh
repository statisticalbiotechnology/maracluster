# change directory to the build directory
cd $1

linux_qt=qtbase-everywhere-src-5.11.2
if [ ! -d Qt-dynamic ]; then

  if [ ! -f ${linux_qt}.tar.xz ]; then
    wget http://download.qt.io/official_releases/qt/5.11/5.11.2/submodules/${linux_qt}.tar.xz
  fi

  tar xf ${linux_qt}.tar.xz
  cd ${linux_qt}

  ./configure -prefix ${build_dir}/tools/Qt-dynamic -opensource -confirm-license -nomake tools -nomake examples -nomake tests
  
  echo "Building Qt base, this may take some time.."
  
  make -j4 > ../qt_installation.log 2>&1
  make install -j4
fi
