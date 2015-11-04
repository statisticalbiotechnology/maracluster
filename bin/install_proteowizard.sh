#!/bin/bash

echo "Download source code for ProteoWizard from their SVN repository"
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
lib_dir=${script_dir}/../libs
#sudo apt-get -y install subversion

mkdir -p ${lib_dir}
svn co -r 7692 https://svn.code.sf.net/p/proteowizard/code/trunk/pwiz ${lib_dir}/proteowizard

# install and keep libraries in the libs folder of this project for linking
cd ${lib_dir}/proteowizard
./clean.sh

echo "Building ProteoWizard and Boost, this may take some time.."

./quickbuild.sh --without-binary-msdata \
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
#cd ../../

