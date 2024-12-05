#!/bin/sh


# CHANGE ME as needed
export CC=/usr/bin/clang
export CXX=/usr/bin/clang++

# Create directories and pull openfhe
rm -rf ./libs/*
[ ! -d "./libs/openfhe" ] && mkdir -p "./libs/openfhe"
[ ! -d "./openfhe" ] && git clone https://github.com/openfheorg/openfhe-development.git openfhe;
cd openfhe || exit;
git checkout 4ebb28e;
# Apply patch, change lib_noiseless.patch to lib.patch as needed
# lib_noiseless sets all standard deviations to 0, ideal for debugging
cp ../patches/* ./src/binfhe;
cd ./src/binfhe || exit;
patch -s -p0 < include.patch;
patch -s -p0 < lib_noiseless.patch;
cd ../../;
# Build openfhe
# Use nativesize = 32 for Parameter sets 1 and 2
[ ! -d "./build" ] && mkdir build;
cd build || exit;
if [  $1 -eq 32 ]; then
  cmake .. -DNATIVE_SIZE=32 -DWITH_NATIVEOPT=ON -DBUILD_STATIC=ON -DBUILD_SHARED=OFF -DWITH_OPENMP=OFF -DCMAKE_INSTALL_PREFIX="$PWD/../../libs/openfhe"
else
  cmake .. -DNATIVE_SIZE=64 -DWITH_NATIVEOPT=ON -DBUILD_STATIC=ON -DBUILD_SHARED=OFF -DWITH_OPENMP=OFF -DCMAKE_INSTALL_PREFIX="$PWD/../../libs/openfhe"
fi
make -j 8;
make install
cd ../../;
