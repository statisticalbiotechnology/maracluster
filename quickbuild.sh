#!/bin/bash
# managing input arguments

src_dir=$(pwd)
build_dir=$src_dir/../build/ubuntu
release_dir=$src_dir/../release/ubuntu

mkdir -p $build_dir/maracluster
#-----cmake-----
cd $build_dir/maracluster;
echo -n "cmake maracluster.....";
cmake -DTARGET_ARCH=amd64 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_PREFIX_PATH=$src_dir/libs $src_dir;
#-----make------
echo -n "make maracluster (this will take few minutes).....";
make -j 4;
make -j 4 package;
sudo make install;

mkdir -p $release_dir
cp -v $build_dir/maracluster/mar*.deb $release_dir
