tools_dir=$1

# change directory to the <build_dir>/tools directory
cd ${tools_dir}

echo "Download source code for ProteoWizard from their TeamCity server"
wget --no-verbose https://teamcity.labkey.org/guestAuth/repository/download/bt81/.lastSuccessful/VERSION
linux_pwiz=pwiz-src-without-tv-$(cat VERSION | sed 's/ /_/g')
# https://teamcity.labkey.org/viewType.html?buildTypeId=bt81
# without-tv: without tests and vendor reader
wget --no-verbose https://teamcity.labkey.org/guestAuth/repository/download/bt81/.lastSuccessful/${linux_pwiz}.tar.bz2

mkdir proteowizard
tar xf ${linux_pwiz}.tar.bz2 --directory proteowizard
cd proteowizard

toolset=""
if [[ "$OSTYPE" == "darwin"* ]]; then
  toolset="toolset=clang"
fi

echo "Building ProteoWizard and Boost, this may take some time.."

# if you have more than 4GB of memory available, you could try to use more than 2 cores to speed things up
./quickbuild.sh ${toolset} -j2 --without-binary-msdata --prefix=../ \
                pwiz/data/common//pwiz_data_common \
                pwiz/data/identdata//pwiz_data_identdata \
                pwiz/data/identdata//pwiz_data_identdata_version \
                pwiz/data/msdata//pwiz_data_msdata_core \
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
                libraries \
                > ../pwiz_installation.log 2>&1

# manually copy some libraries and headers used by maracluster but not by proteowizard
find build-*-x86_64/ -type f | grep -i libboost_regex-.*\.a$ | xargs -I{} cp {} ../lib
find build-*-x86_64/ -type f | grep -i libboost_program_options-.*\.a$ | xargs -I{} cp {} ../lib

rsync -ap --include "*/" --include "*.h" --include "*.hpp" --exclude "*"  libraries/zlib-1.2.3/ ../include
rsync -ap --include "*/" --include "*.h" --include "*.hpp" --exclude "*" libraries/boost_aux/boost/ ../include/boost
rsync -ap --include "*/" --include '*.ipp' --exclude '*' libraries/boost_1_76_0/boost/ ../include/boost

cat ../pwiz_installation.log
ls ../lib

# the boost libraries' naming convention does not always work well with cmake, so we force a more simple naming convention
for component in system thread chrono regex filesystem iostreams program_options serialization; do
  ln -s -f $(ls ../lib/libboost_${component}-*.a | head -n1) ../lib/libboost_${component}.a
done

# copy the boost::asio library, which is not included by the ProteoWizard boost tar but is needed for maracluster
cd ${tools_dir}
wget --no-verbose --no-check-certificate https://sourceforge.net/projects/asio/files/asio/1.18.2%20%28Stable%29/boost_asio_1_18_2.tar.gz/download -O boost_asio.tar.gz
tar -xzf boost_asio.tar.gz
rsync -ap --include "*/" --include "*.hpp" --include '*.ipp' --exclude '*' boost_asio_1_18_2/boost/ include/boost
