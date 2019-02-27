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
    other_packages="cmake gnutar wget coreutils"
elif [[ -f ${HOME}/bin/brew ]]
  then
    echo "[ Package manager ] : Homebrew "
    package_manager=$HOME/bin/brew
    other_packages="cmake gnu-tar wget coreutils"
elif [[ -f /usr/local/bin/brew ]]
  then
    echo "[ Package manager ] : Homebrew "
    package_manager="brew"
    ${package_manager} update || true # brew.rb raises an error on the vagrant box, just ignore it
    other_packages="cmake gnu-tar wget coreutils"
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

rm -f $build_dir/{maracluster,maracluster-gui}/mar*.dmg

echo "The Builder $0 is building the MaRaCluster packages with src=${src_dir} an\
d build=${build_dir} for user" `whoami`
$package_manager install $other_packages
CMAKE_BINARY=cmake # this can be overridden if a newer version of cmake is needed

#----------------------------------------

mkdir -p ${build_dir}/tools
cd ${build_dir}/tools

if [ ! -d ${build_dir}/tools/proteowizard ]; then
  ${src_dir}/maracluster/admin/builders/install_proteowizard.sh ${build_dir}/tools
fi

#-----MaRaCluster-GUI dependencies-------

if [ "$no_gui" != true ] ; then
  source ${src_dir}/maracluster/admin/builders/install_qt.sh ${build_dir}/tools
fi

#-------------------------------------------

mkdir -p $build_dir/maracluster
#-----cmake-----
# we need to install to /usr/local instead of /usr: https://github.com/Benjamin-Dobell/Heimdall/issues/291
cd $build_dir/maracluster;
echo -n "cmake maracluster.....";
${CMAKE_BINARY} -DCMAKE_CXX_COMPILER="/usr/bin/clang++" -DTARGET_ARCH="x86_64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_PREFIX_PATH="${build_dir}/tools/" $src_dir/maracluster;

#-----make------
echo -n "make maracluster.....";
make -j 4;
make -j 4 package;

echo "build directory was : ${build_dir}";

mkdir -p $release_dir
cp -v $build_dir/maracluster/mar*.dmg $release_dir

if [ "$no_gui" != true ] ; then
  #######maracluster-gui########
  mkdir -p $build_dir/maracluster-gui
  cd $build_dir/maracluster-gui
  
  ${CMAKE_BINARY} -DCMAKE_CXX_COMPILER="/usr/bin/clang++" -DTARGET_ARCH="x86_64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_PREFIX_PATH="${build_dir}/tools/;${build_dir}/tools/Qt-dynamic" $src_dir/maracluster/src/qt-gui/
  
  echo -n "make maracluster-gui.....";
  make -j4
  make -j4 package;
  
  cp -v $build_dir/maracluster-gui/mar*.dmg $release_dir
fi

#--------------------------------------------



