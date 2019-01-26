#!/bin/bash
# Requirements are:
# XCode
# Command line tools (check if installed with "xcode-select -p", otherwise install with "xcode-select --install")
# MacPorts or homebrew as package manager
# PackageMaker (https://developer.apple.com/downloads ::search::
#----------------------------------------

# managing input arguments
while getopts “s:b:r:t:” OPTION; do
  case $OPTION in
    s) src_dir=${OPTARG};;
    t) branch=${OPTARG};;
    r) release_dir=${OPTARG};;
    b) build_dir=${OPTARG};;
    g) no_gui=true;;
    \?) echo "Invalid option: -${OPTARG}" >&2;;
  esac
done

if [[ ! -d /Applications/XCode.app ]]
  then
    echo "Apple developer tools are required (Search for XCode in the App Store)"
    exit 1
fi

if [[ ! -d /Applications/PackageMaker.app ]]
  then
    echo "Apple developer PackageManager is required and expected in the "
    echo "/Applications folder. If you have moved it elsewhere, please change this script"
    echo ""
    echo "It is part of the Auxiliary tools for XCode - Late July 2012"
    echo "Yes, 2012! since then Apple moved to the app store and requires"
    echo "packages and dmgs to be build differently. "
    echo "However, the old packagemaker still works with 10.11"
    echo
    echo "You can find it here: "
    echo "http://adcdownload.apple.com/Developer_Tools/auxiliary_tools_for_xcode__late_july_2012/xcode44auxtools6938114a.dmg"
    echo ""
    exit 1
fi

package_manager_installed=true
if [[ -d /opt/local/var/macports ]]
  then
    echo "[ Package manager ] : MacPorts "
    package_manager="sudo port"
    other_packages="cmake gnutar"
elif [[ -f ${HOME}/bin/brew ]]
  then
    echo "[ Package manager ] : Homebrew "
    package_manager=$HOME/bin/brew
    other_packages="cmake gnutar"
else
    package_manager_installed=false
fi

if [ "$package_manager_installed" == false ]
  then
  echo "Error: no suitable package manager installed"
  echo " Get homebrew or macports:"
  echo "  Homebrew: http://brew.sh/ "
  echo "  MacPorts: http://www.macports.org/install.php"
  exit 1
fi

if [[ -z ${build_dir} ]]; then
  build_dir="$(mktemp -d -t build)";
fi
if [[ -z ${src_dir} ]]; then
  if [[ -n  ${branch} ]]; then
    if [[ ! -f /usr/bin/git ]]; then
      $package_manager install git;
    fi
    src_dir="$(mktemp -d -t src)";
    git clone --branch "$1" https://github.com/statisticalbiotechnology/maracluster.git "${src_dir}/maracluster";
  else
    # Might not work if we have symlinks in the way
    src_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )/../../../" && pwd )
  fi
fi
if [[ -z ${release_dir} ]]; then
  release_dir=${HOME}/release
fi

rm $build_dir/maracluster/mar*.dmg

echo "The Builder $0 is building the MaRaCluster packages with src=${src_dir} an\
d build=${build_dir} for user" `whoami`
$package_manager install $other_packages

#----------------------------------------

mkdir -p ${build_dir}/tools
cd ${build_dir}/tools

if [ ! -d ${build_dir}/tools/proteowizard ]; then
  echo "Download source code for ProteoWizard from their SVN repository"
  $package_manager install subversion
  rev=9393
  svn co -r ${rev} --depth immediates svn://svn.code.sf.net/p/proteowizard/code/trunk/pwiz ./proteowizard
  svn update -r ${rev} --set-depth infinity ./proteowizard/pwiz
  svn update -r ${rev} --set-depth infinity ./proteowizard/libraries

  # install and keep libraries in the libs folder of this project for linking
  cd proteowizard

  ./clean.sh

  echo "Building ProteoWizard and Boost, this may take some time.."
  
  # if you have more than 4GB of memory available, you could try to use more than 2 cores to speed things up
  ./quickbuild.sh -j4 --without-binary-msdata \
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
  find build-macosx-x86_64/ -type f | grep -i .a$ | xargs -I{} cp {} ../lib
  
  mkdir -p ../include
  rsync -ap --include "*/" --include "*.h" --include "*.hpp" --include "*.ipp" --exclude "*" pwiz libraries/boost_1_56_0/boost libraries/boost_aux/boost ../include/
  rsync -ap --include "*/" --include "*.h" --include "*.hpp" --exclude "*"  libraries/zlib-1.2.3/ ../include/zlib
fi

#-----MaRaCluster-GUI dependencies-------

if [ "$no_gui" != true ] ; then
  if [ ! -d ${build_dir}/tools/Qt-dynamic ]; then
    cd ${build_dir}/tools

    if [ ! -f qtbase-everywhere-src-5.11.2.tar.xz ]; then
      wget http://download.qt.io/official_releases/qt/5.11/5.11.2/submodules/qtbase-everywhere-src-5.11.2.tar.xz
    fi

    tar xf qtbase-everywhere-src-5.11.2.tar.xz
    cd qtbase-everywhere-src-5.11.2

    ./configure -prefix ../build/Qt-dynamic -opensource -confirm-license -nomake tools -nomake examples -nomake tests

    make -j4
    make install -j4
  fi
fi

#-------------------------------------------

mkdir -p $build_dir/maracluster
#-----cmake-----
# we need to install to /usr/local instead of /usr: https://github.com/Benjamin-Dobell/Heimdall/issues/291
cd $build_dir/maracluster;
echo -n "cmake maracluster.....";
cmake -DCMAKE_CXX_COMPILER="/usr/bin/clang++" -DTARGET_ARCH="x86_64" -DBOOST_ROOT="${build_dir}/tools/proteowizard/libraries/boost_1_56_0" -DCMAKE_BUILD_TYPE=Release -DBoost_COMPILER=-xgcc42 -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_PREFIX_PATH="${build_dir}/tools/" $src_dir/maracluster;
#-----make------
echo -n "make maracluster (this will take few minutes).....";
make -j 4;
make -j 4 package;
#sudo make install;

echo "build directory was : ${build_dir}";

mkdir -p $release_dir
cp -v $build_dir/maracluster/mar*.dmg $release_dir

if [ "$no_gui" != true ] ; then
  #######maracluster-gui########
  mkdir -p $build_dir/maracluster-gui
  cd $build_dir/maracluster-gui
  
  cmake -DCMAKE_CXX_COMPILER="/usr/bin/clang++" -DTARGET_ARCH="x86_64" -DBOOST_ROOT="${build_dir}/tools/proteowizard/libraries/boost_1_56_0" -DCMAKE_BUILD_TYPE=Release -DBoost_COMPILER=-xgcc42 -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_PREFIX_PATH="${build_dir}/tools/;${build_dir}/tools/Qt-dynamic" $src_dir/maracluster/src/qt-gui/
  make -j4
  make -j4 package
  
  cp -v $build_dir/maracluster-gui/mar*.dmg $release_dir
fi

#--------------------------------------------



