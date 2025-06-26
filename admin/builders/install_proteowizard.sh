tools_dir=$1

# change directory to the <build_dir>/tools directory
cd ${tools_dir}

echo "Download source code for ProteoWizard from their TeamCity server"
wget --no-check-certificate --no-verbose -O bt81.xml https://proteowizard.sourceforge.io/releases/bt81.xml
# without-tv: without tests and vendor reader
# the bt81.xml is formatted to contain all content on a single line, so we use a regex that matches the filename until .tar.bz2 but before the closing xml tag ("[^<]*")
read BUILD_ID FILE_NAME < <(sed -n 's/.*id:\([0-9]*\)\/artifacts\/content\/\(pwiz-src-without-tv-[^<]*\.tar\.bz2\).*/\1 \2/p' bt81.xml)

if [ ! -f ${tools_dir}/${FILE_NAME} ]; then
  wget --no-check-certificate --no-verbose -O ${FILE_NAME} https://mc-tca-01.s3.us-west-2.amazonaws.com/ProteoWizard/bt81/${BUILD_ID}/${FILE_NAME}

  mkdir proteowizard
  tar xf ${FILE_NAME} --directory proteowizard  
fi

cd proteowizard

# undo position independent code fix introduced here: https://github.com/ProteoWizard/pwiz/pull/1980
sed -i.bak 's/ -no-pie -fno-pie//g' Jamroot.jam

echo "Building ProteoWizard and Boost, this may take some time.."

# if you have more than 4GB of memory available, you could try to use more than 2 cores to speed things up
# add -d2 flag to get verbose output with compiler and linker commands
./quickbuild.sh ${toolset} -j2 --prefix=../ \
                pwiz/data/common//pwiz_data_common \
                pwiz/data/identdata//pwiz_data_identdata \
                pwiz/data/identdata//pwiz_data_identdata_version \
                pwiz/data/msdata//pwiz_data_msdata_core \
                pwiz/data/msdata//pwiz_data_msdata \
                pwiz/data/msdata//pwiz_data_msdata_version \
                pwiz/data/msdata/mzmlb//pwiz_data_msdata_mzmlb \
                pwiz/data/proteome//pwiz_data_proteome \
                pwiz/utility/chemistry//pwiz_utility_chemistry \
                pwiz/utility/minimxml//pwiz_utility_minimxml \
                pwiz/utility/misc//SHA1 \
                pwiz/utility/misc//pwiz_utility_misc \
                /ext/zlib//z \
                /ext/hdf5//hdf5 \
                /ext/hdf5//hdf5pp \
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

status=$?
if [ $status -ne 0 ]; then
    echo "❌ Build failed with status $status. Showing log:"
    cat ../pwiz_installation.log
    exit $status
fi

# manually copy some libraries and headers used by maracluster but not by proteowizard
find build-*-x86_64/ -type f | grep -i libboost_regex-.*\.a$ | xargs -I{} cp {} ../lib
find build-*-x86_64/ -type f | grep -i libboost_program_options-.*\.a$ | xargs -I{} cp {} ../lib

rsync -ap --include "*/" --include "*.h" --include "*.hpp" --exclude "*"  libraries/zlib-1.2.3/ ../include
rsync -ap --include "*/" --include "*.h" --include "*.hpp" --exclude "*" libraries/boost_aux/boost/ ../include/boost
rsync -ap --include "*/" --include '*.ipp' --exclude '*' libraries/boost_1_86_0/boost/ ../include/boost

ls ../lib

# the boost libraries' naming convention does not always work well with cmake, so we force a more simple naming convention
for component in system thread chrono regex filesystem iostreams program_options serialization; do
  ln -s -f $(ls ../lib/libboost_${component}-*.a | head -n1) ../lib/libboost_${component}.a
done

# copy the boost::asio and boost::unordered library, which are not included by the ProteoWizard boost tar but are needed for maracluster
cd ${tools_dir}
wget --no-verbose --no-check-certificate https://sourceforge.net/projects/asio/files/asio/1.18.2%20%28Stable%29/boost_asio_1_18_2.tar.gz/download -O boost_asio.tar.gz
tar -xzf boost_asio.tar.gz
rsync -ap --include "*/" --include "*.hpp" --include '*.ipp' --exclude '*' boost_asio_1_18_2/boost/ include/boost

wget --no-verbose --no-check-certificate https://github.com/boostorg/unordered/archive/refs/tags/boost-1.86.0.tar.gz -O boost_unordered.tar.gz
tar -xzf boost_unordered.tar.gz
rsync -ap --include "*/" --include "*.hpp" --include '*.ipp' --exclude '*' unordered-boost-1.86.0/include/boost/ include/boost

touch pwiz_successful.txt