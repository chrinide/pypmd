#!/bin/bash

GCC5ROOT=/home/jluis/pkgs/gcc/6.3.0-cuda
CUDA=/usr/local/cuda

mkdir -p $GCC5ROOT/source
cd $GCC5ROOT/source
git clone https://github.com/MentorEmbedded/nvptx-newlib.git
git clone https://github.com/MentorEmbedded/nvptx-tools.git
wget http://gcc.parentingamerica.com/releases/gcc-6.3.0/gcc-6.3.0.tar.bz2 
tar jxvf gcc-6.3.0.tar.bz2
mv gcc-6.3.0 gcc

mkdir -p $GCC5ROOT/build/nvptx-build
cd $GCC5ROOT/build/nvptx-build
$GCC5ROOT/source/nvptx-tools/configure  \
  --prefix=$GCC5ROOT \
  --target=nvptx-none \
  --with-cuda-driver-include=$CUDA/include \
  --with-cuda-driver-lib=/usr/lib64/nvidia \
  --with-cuda-runtime-include=$CUDA/include \
  --with-cuda-runtime-lib=$CUDA/lib64 \
  CC='gcc -m64' \
  CXX='g++ -m64'
make -j12
make install

mkdir -p $GCC5ROOT/build/gcc6-accel
cd $GCC5ROOT/build/gcc6-accel
ln -vs $GCC5ROOT/source/nvptx-newlib/newlib $GCC5ROOT/source/gcc/newlib
ln -vs . $GCC5ROOT/install/nvptx-none/usr
target=$($GCC5ROOT/source/gcc/config.guess)
$GCC5ROOT/source/gcc/configure \
  --prefix= \
  --target=nvptx-none \
  --enable-as-accelerator-for="$target" \
  --enable-languages=c,c++,fortran,lto \
  --enable-checking=yes,df,fold,rtl \
  --disable-multilib \
  --with-sysroot=/nvptx-none \
  --with-build-sysroot=$GCC5ROOT/nvptx-none \
  --with-build-time-tools=$GCC5ROOT/nvptx-none/bin \
  --disable-sjlj-exceptions \
  --enable-newlib-io-long-long \
  CC='gcc -m64'\
  CXX='g++ -m64'
make -j12
make DESTDIR=$GCC5ROOT install

mkdir -p $GCC5ROOT/build/gcc6
cd $GCC5ROOT/build/gcc6
$GCC5ROOT/source/gcc/configure \
  --prefix= \
  --disable-bootstrap \
  --enable-languages=c,c++,fortran,lto \
  --disable-multilib \
  --enable-offload-targets=nvptx-none=$GCC5ROOT \
  --with-cuda-driver-include=$CUDA/include \
  CC='gcc -m64' \
  CXX='g++ -m64' \
  --with-sysroot=

make -j12
make DESTDIR=$GCC5ROOT install
