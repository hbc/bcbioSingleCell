#!/bin/bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    # homebrew bottle
    brew update
    brew install hdf5
else
    # build from source
    wd=$PWD
    cd ..
    wget "$HDF5_RELEASE_URL/hdf5-${HDF5_VERSION%.*}/hdf5-$HDF5_VERSION/src/hdf5-$HDF5_VERSION.tar.gz"
    tar -xzf "hdf5-$HDF5_VERSION.tar.gz"
    cd "hdf5-$HDF5_VERSION"
    ./configure --prefix=/usr/local
    sudo make install
    cd $wd
fi
