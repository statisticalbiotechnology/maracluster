#!/bin/bash
# managing input arguments
while getopts “s:b:r:t:g” OPTION; do
  case $OPTION in
    s) src_dir=${OPTARG};;
    t) branch=${OPTARG};;
    r) release_dir=${OPTARG};;
    b) build_dir=${OPTARG};;
    g) no_gui=true;;
    \?) echo "Invalid option: -${OPTARG}" >&2;;
  esac
done

if [[ -z ${build_dir} ]]; then
  build_dir="$(mktemp -d --tmpdir build_XXXX)";
fi
if [[ -z ${src_dir} ]]; then
  if [[ -n  ${branch} ]]; then
    sudo apt-get install git;
    src_dir="$(mktemp -d --tmpdir build_XXXX)";
    git clone --branch "$1" https://github.com/statisticalbiotechnology/maracluster.git "${src_dir}/maracluster";
  else
    src_dir=$(dirname ${BASH_SOURCE})/../../../
  fi
fi
if [[ -z ${release_dir} ]]; then
  release_dir=${HOME}/release
fi

rm -f $build_dir/{maracluster,maracluster-gui}/mar*.rpm

sudo dnf install -y gcc gcc-c++ rpm-build cmake rsync
CMAKE_BINARY=cmake # this can be overridden if a newer version of cmake is needed

mkdir -p ${build_dir}/tools
cd ${build_dir}/tools

# we need gcc >= 5 to compile ProteoWizard, since they use the c++14 flag
pre=""
if [[ $(rpm -q --queryformat '%{VERSION}' centos-release) < 8 ]]; then
  sudo dnf install -y centos-release-scl
  sudo dnf install -y devtoolset-7-gcc*
  pre="scl enable devtoolset-7 --"
fi

if [ ! -f ${build_dir}/tools/pwiz_successful.txt ]; then
  ${pre} ${src_dir}/maracluster/admin/builders/install_proteowizard.sh ${build_dir}/tools
fi

#-----MaRaCluster-GUI dependencies-------

if [ "$no_gui" != true ] ; then
  # patchelf is available from Ubuntu 16.04, install from source otherwise
  sudo dnf install -y patchelf || missing_patchelf=true
  if [ "$missing_patchelf" == true ] ; then
    if [ ! -d ${build_dir}/tools/patchelf-0.8 ]; then
      wget --no-check-certificate http://nixos.org/releases/patchelf/patchelf-0.8/patchelf-0.8.tar.bz2
      tar xf patchelf-0.8.tar.bz2
      cd patchelf-0.8/
      ./configure --prefix="${build_dir}/tools"
      make install
    fi
  fi
  
  sudo dnf install -y mesa-libGL-devel libicu-devel freetype-devel
  ${pre} source ${src_dir}/maracluster/admin/builders/install_qt.sh ${build_dir}/tools
  
  # the child process in devtoolset-7 does not export the CMAKE_BINARY variable from install_qt.sh, so set it again to be safe
  function version_lt() { test "$(echo "$@" | tr " " "\n" | (sort -rV || gsort -rV) | head -n 1)" != "$1"; }

  if version_lt $(cmake --version | head -n1 | cut -f3 -d ' ') "3.5"; then
    CMAKE_BINARY=${build_dir}/tools/bin/cmake
  fi
fi

mkdir -p $build_dir/maracluster
#-----cmake-----
cd $build_dir/maracluster;
echo -n "cmake maracluster.....";
${pre} ${CMAKE_BINARY} -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH=$build_dir/tools $src_dir/maracluster;
#-----make------
echo -n "make maracluster (this will take few minutes).....";
${pre} make -j 4;
${pre} make -j 4 package;

if [ "$no_gui" != true ] ; then
  #######maracluster-gui########
  mkdir -p $build_dir/maracluster-gui
  cd $build_dir/maracluster-gui
  #-----cmake-----
  echo -n "cmake maracluster-gui.....";
  ${pre} ${CMAKE_BINARY} -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH="$build_dir/tools/;$build_dir/tools/Qt-dynamic/" $src_dir/maracluster/src/qt-gui;

  #-----make------
  echo -n "make maracluster-gui (this will take few minutes).....";

  ${pre} make -j 4;
  ${pre} make -j 4 package;
fi

mkdir -p $release_dir
cp -v $build_dir/maracluster/mar*.rpm $release_dir && \
  ([ "$no_gui" == true ] || cp -v $build_dir/maracluster-gui/mar*.rpm $release_dir)


