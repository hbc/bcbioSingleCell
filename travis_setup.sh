#!/bin/bash

if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    wd=$PWD
    cd ..
    wget "$HDF5_RELEASE_URL/hdf5-${HDF5_VERSION%.*}/hdf5-$HDF5_VERSION/src/hdf5-$HDF5_VERSION.tar.gz"
    tar -xzf "hdf5-$HDF5_VERSION.tar.gz"
    cd "hdf5-$HDF5_VERSION"
    ./configure --prefix=/usr/local
    sudo make install
    cd $wd
fi
