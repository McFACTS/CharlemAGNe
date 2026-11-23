#!/bin/bash

SOURCE_DIR="$(dirname "$(readlink -f "$0")")/.."

### Build HDF5 ###
# If we're building HDF5, do it after we set *all flags*

H5VER=1.14.2
H5VERU=1_14_2

# Download, blow away existing dir on option
cd $SOURCE_DIR/external
# Allow complete reconfigure (for switching compilers, takes longer)
# TODO allow not
rm -rf hdf5-${H5VER}/

# Download if needed
if [ ! -f hdf5-${H5VER}.tar.gz ]; then
  curl https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_${H5VERU}/src/hdf5-${H5VER}.tar.gz -o hdf5-${H5VER}.tar.gz
fi
# Unpack if needed (or deleted)
if [ ! -d hdf5-${H5VER}/ ]; then
  tar xf hdf5-${H5VER}.tar.gz
fi

# Configure and compile
cd $SOURCE_DIR/external/hdf5-${H5VER}/

echo Configuring HDF5...

export CFLAGS="-fPIC $CFLAGS"
sh ./configure -C --prefix=$SOURCE_DIR/external/hdf5 --enable-build-mode=production \
  --disable-dependency-tracking --disable-tests --disable-tools --disable-deprecated-symbols \
  --disable-shared --without-zlib > build-hdf5.log

sleep 1

echo "Building HDF5 (probably 30s-2min)"
# Should be fine anywhere /shrug
make -j >> build-hdf5.log 2>&1
make install >> build-hdf5.log 2>&1
make clean >> build-hdf5.log 2>&1
cd $SOURCE_DIR

echo Built HDF5 version $H5VER
