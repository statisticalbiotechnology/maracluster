# change directory to the build directory
cd $1

echo "Download source code for ProteoWizard from their TeamCity server"
linux_pwiz=pwiz-src-without-tv-3_0_19025_7f0e41d
# https://teamcity.labkey.org/viewType.html?buildTypeId=bt81
# without-tv: without tests and vendor readers
wget https://teamcity.labkey.org/guestAuth/repository/download/bt81/.lastSuccessful/${linux_pwiz}.tar.bz2

mkdir proteowizard
tar xf ${linux_pwiz}.tar.bz2 --directory proteowizard
cd proteowizard

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

# install and keep libraries in the libs folder of this project for linking
# for revision 7692 the "libraries" target does not work, so we have to copy everything manually
mkdir -p ../lib
find build-linux-x86_64/ -type f | grep -i .a$ | xargs -i cp {} ../lib

# the boost libraries' naming convention does not always work well with cmake, so we force a more simple naming convention
ln -s -f ../lib/libboost_system-*.a ../lib/libboost_system.a
ln -s -f ../lib/libboost_thread-*.a ../lib/libboost_thread.a
ln -s -f ../lib/libboost_chrono-*.a ../lib/libboost_chrono.a
ln -s -f ../lib/libboost_regex-*.a ../lib/libboost_regex.a
ln -s -f ../lib/libboost_filesystem-*.a ../lib/libboost_filesystem.a
ln -s -f ../lib/libboost_iostreams-*.a ../lib/libboost_iostreams.a
ln -s -f ../lib/libboost_program_options-*.a ../lib/libboost_program_options.a
ln -s -f ../lib/libboost_serialization-*.a ../lib/libboost_serialization.a

mkdir -p ../include
find pwiz/ -type f | grep -i '.h$\|.hpp$' | xargs -i cp --parents {} ../include/

cd libraries/boost_1_56_0/
find boost/ -type f | grep -i '.h$\|.hpp$' | xargs -i cp --parents {} ../../../include/

cd ../zlib-1.2.3
find ./ -type f | grep -i '.h$\|.hpp$' | xargs -i cp --parents {} ../../../include/
