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

#rm -f $build_dir/maracluster/mar*.deb

sudo apt-get update;
sudo apt-get upgrade;
sudo apt-get -y install g++ make cmake

mkdir -p ${build_dir}/tools
cd ${build_dir}/tools

if [ ! -d ${build_dir}/tools/proteowizard ]; then
  echo "Download source code for ProteoWizard from their SVN repository"
  sudo apt-get -y install subversion
  rev=9393
  svn co -r ${rev} --depth immediates svn://svn.code.sf.net/p/proteowizard/code/trunk/pwiz ./proteowizard
  svn update -r ${rev} --set-depth infinity ./proteowizard/pwiz
  svn update -r ${rev} --set-depth infinity ./proteowizard/libraries

  # install and keep libraries in the libs folder of this project for linking
  cd proteowizard

  ./clean.sh

  echo "Building ProteoWizard and Boost, this may take some time.."
  
  # if you have more than 4GB of memory available, you could try to use more than 2 cores to speed things up
  ./quickbuild.sh -j2 --without-binary-msdata \
                  pwiz/data/common//pwiz_data_common \
                  pwiz/data/identdata//pwiz_data_identdata \
                  pwiz/data/identdata//pwiz_data_identdata_version \
                  pwiz/data/msdata//pwiz_data_msdata \
                  pwiz/data/msdata//pwiz_data_msdata_version \
                  pwiz/data/proteome//pwiz_data_proteome \
                  pwiz/utility/chemistry//pwiz_utility_chemistry \
                  pwiz/utility/minimxml//pwiz_utility_minimxml \
                  pwiz/utility/misc//SHA1 \
                  pwiz/utility/misc//pwiz_utility_misc \
                  /ext/zlib//z \
                  /ext/boost//system \
                  /ext/boost//thread \
                  /ext/boost//chrono \
                  /ext/boost//regex \
                  /ext/boost//filesystem \
                  /ext/boost//iostreams \
                  /ext/boost//program_options \
                  /ext/boost//serialization \
                  > ../pwiz_installation.log 2>&1

  # for revision 7692 the "libraries" target does not work, so we have to copy everything manually
  mkdir -p ../lib
  find build-linux-x86_64/ -type f | grep -i .a$ | xargs -i cp {} ../lib
  
  mkdir -p ../include
  find pwiz/ -type f | grep -i '.h$\|.hpp$' | xargs -i cp --parents {} ../include/
  
  cd libraries/boost_1_56_0/
  find boost/ -type f | grep -i '.h$\|.hpp$' | xargs -i cp --parents {} ../../../include/
  
  cd ../zlib-1.2.3
  find ./ -type f | grep -i '.h$\|.hpp$' | xargs -i cp --parents {} ../../../include/
fi

#-----MaRaCluster-GUI dependencies-------

if [ "$no_gui" != true ] ; then
  # patchelf is available from Ubuntu 16.04, install from source otherwise
  sudo apt -y install patchelf || missing_patchelf=true
  if [ "$missing_patchelf" == true ] ; then
    wget http://nixos.org/releases/patchelf/patchelf-0.8/patchelf-0.8.tar.bz2
    tar xf patchelf-0.8.tar.bz2
    cd patchelf-0.8/
    ./configure --prefix="${build_dir}/tools"
    make install
    patchelf_binary=${build_dir}/tools/bin/patchelf
  else
    patchelf_binary=patchelf
  fi
  
  ubuntu_qt=qtbase-everywhere-src-5.11.2
  if [ ! -d ${build_dir}/tools/Qt-dynamic ]; then
    cd ${build_dir}/tools

    if [ ! -f ${ubuntu_qt}.tar.xz ]; then
      wget http://download.qt.io/official_releases/qt/5.11/5.11.2/submodules/${ubuntu_qt}.tar.xz
    fi

    tar xf ${ubuntu_qt}.tar.xz
    cd ${ubuntu_qt}

    ./configure -prefix ${build_dir}/tools/Qt-dynamic -opensource -confirm-license -nomake tools -nomake examples -nomake tests

    make -j4
    make install -j4
  fi
fi

mkdir -p $build_dir/maracluster
#-----cmake-----
cd $build_dir/maracluster;
echo -n "cmake maracluster.....";
cmake -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH=$build_dir/tools $src_dir/maracluster;
#-----make------
echo -n "make maracluster (this will take few minutes).....";
make -j 4;
make -j 4 package;
sudo make install;

mkdir -p $release_dir
cp -v $build_dir/maracluster/mar*.deb $release_dir

if [ "$no_gui" != true ] ; then
  #######maracluster-gui########
  mkdir -p $build_dir/maracluster-gui
  cd $build_dir/maracluster-gui
  #-----cmake-----
  echo -n "cmake maracluster-gui.....";
  cmake -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH="$build_dir/tools/;$build_dir/tools/Qt-dynamic/" $src_dir/maracluster/src/qt-gui;

  #-----make------
  echo -n "make maracluster-gui (this will take few minutes).....";

  make -j 4;
  make -j 4 package;
  sudo make install
  
  cp -v $build_dir/maracluster-gui/mar*.deb $release_dir
fi

