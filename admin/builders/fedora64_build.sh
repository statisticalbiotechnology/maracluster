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

sudo dnf install -y gcc gcc-c++ rpm-build cmake
CMAKE_BINARY=cmake # this can be overridden if a newer version of cmake is needed

mkdir -p ${build_dir}/tools
cd ${build_dir}/tools

if [ ! -f ${build_dir}/tools/pwiz_successful.txt ]; then
  ${src_dir}/maracluster/admin/builders/install_proteowizard.sh ${build_dir}/tools
fi

#-----MaRaCluster-GUI dependencies-------

if [ "$no_gui" != true ] ; then
  sudo dnf install -y patchelf mesa-libGL-devel libicu-devel freetype-devel
  source ${src_dir}/maracluster/admin/builders/install_qt.sh ${build_dir}/tools
fi

mkdir -p $build_dir/maracluster
#-----cmake-----
cd $build_dir/maracluster;
echo -n "cmake maracluster.....";
${CMAKE_BINARY} -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH=$build_dir/tools $src_dir/maracluster;
#-----make------
echo -n "make maracluster (this will take few minutes).....";
make -j 4;
make -j 4 package;

if [ "$no_gui" != true ] ; then
  #######maracluster-gui########
  mkdir -p $build_dir/maracluster-gui
  cd $build_dir/maracluster-gui
  #-----cmake-----
  echo -n "cmake maracluster-gui.....";
  ${CMAKE_BINARY} -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH="$build_dir/tools/;$build_dir/tools/Qt-dynamic/" $src_dir/maracluster/src/qt-gui;

  #-----make------
  echo -n "make maracluster-gui (this will take few minutes).....";

  make -j 4;
  make -j 4 package;
fi

mkdir -p $release_dir
cp -v $build_dir/maracluster/mar*.rpm $release_dir && \
  ([ "$no_gui" == true ] || cp -v $build_dir/maracluster-gui/mar*.rpm $release_dir)

