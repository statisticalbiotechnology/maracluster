#!/bin/bash
# Requirements are:
# XCode
# Command line tools (check if installed with "xcode-select -p", otherwise install with "xcode-select --install")
# MacPorts or homebrew as package manager
# PackageMaker (https://developer.apple.com/downloads ::search::
#----------------------------------------

# managing input arguments
no_gui=false
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

if [[ ! -d /Applications/XCode.app ]]
  then
    echo "Apple developer tools are required (Search for XCode in the App Store)"
    exit 1
fi

package_manager_installed=false
if [[ -d /opt/local/var/macports ]]; then
    echo "[ Package manager ] : MacPorts"
    package_manager="sudo port"
    other_packages="gnutar wget coreutils libomp cmake"
    package_manager_installed=true
else
    brew_paths=(
        "$HOME/bin/brew"
        "/opt/homebrew/bin/brew"
        "/usr/local/bin/brew"
    )

    for path in "${brew_paths[@]}"; do
        if [[ -x "$path" ]]; then
            echo "[ Package manager ] : Homebrew"
            package_manager="$path"
            other_packages="gnu-tar wget coreutils libomp cmake"
            "$package_manager" update || true
            package_manager_installed=true
            break
        fi
    done
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

if [ ! -f ${build_dir}/tools/pwiz_successful.txt ]; then
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

if [ "$no_gui" != true ] ; then
  #######maracluster-gui########
  mkdir -p $build_dir/maracluster-gui
  cd $build_dir/maracluster-gui
  
  ${CMAKE_BINARY} -DCMAKE_CXX_COMPILER="/usr/bin/clang++" -DTARGET_ARCH="x86_64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_PREFIX_PATH="${build_dir}/tools/;${build_dir}/tools/Qt-dynamic" $src_dir/maracluster/src/qt-gui/
  
  echo -n "make maracluster-gui.....";
  make -j4
  make -j4 package;
fi

mkdir -p $release_dir
cp -v $build_dir/maracluster/mar*.pkg $release_dir && \
  ([ "$no_gui" == true ] || cp -v $build_dir/maracluster-gui/mar*.dmg $release_dir)

