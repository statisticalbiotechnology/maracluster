#!/bin/bash
# managing input arguments

src_dir=$(pwd)
build_dir=$src_dir/../build/ubuntu
release_dir=$src_dir/../release/ubuntu

cd admin/builders
./ubuntu64_build.sh -b ${build_dir} -r ${release_dir}
