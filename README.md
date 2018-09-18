# DilateAF: Dilate using ArrayFire and TileDB

This program dilates a large 3D voxel data (representing a foam skeleton structure) that might not fit in memory. It uses TileDB to store voxel data in compressed format on disk and efficiently retreive chunks of data as and when required. ArrayFire is used to perform dilation operation using FFT convolution. The distributive property of morphological dilation is exploited to dilate chunks of structure at a time and then merging (union) the dilated outputs.

### Usage instructions:
```sh
$ cd /path/to/project-directory
$ mkdir build
$ cd build
$ cmake -DArrayFire_DIR=/path/to/arrayfire/build ..
$ make
$ ./dilateAF
```

### Dependencies:
* [ArrayFire] (http://arrayfire.org/)
* [TileDB] (https://tiledb.io/)

#### Installation instructions for ArrayFire on Ubuntu:
* use [ArrayFire (AF) Installation Instructions for Linux](https://github.com/arrayfire/arrayfire/wiki/Build-Instructions-for-Linux) with CPU backend - install all dependencies
* cmake >3.5 is required for AF; install cmake 3.5 or above from source if you don't have it; configure, make, and install
* configure AF using cmake
* `make -j8` (*number of threads*)
* if curl gives error for unsupported protocol, these are few things that you could try:
  - install libcurl from source and configure with `--with-ssl` option, make, and install
  - change the single quotes to double quotes for url in `boost_compute-urlinfo.txt` file
  - download the zip file at the url and rename it appropriately
* run tests and examples

#### Installation instructions for TileDB on Ubuntu:
* use [TileDB Installation Instructions for Linux - Building from Source](https://docs.tiledb.io/en/latest/installation.html)
* install dependencies: 
  -  zlib (zlib1g-dev)
  -  LZ4 (liblz4-dev)
  -  bzip2 (libbz2-dev)
  -  Zstandard (from source)
  -  Blosc (from source)
* git clone and create build dir
* configure using cmake
* build and install using make
* if curl gives error for unsupported protocol, these are few things that you could try:
  - install libcurl from source and configure with `--with-ssl` option, make, and install
  - change single quotes to double quotes for urls of tbb, catch, and spdlog in `./externals/src/....`
  - download the zip files from the urls and rename them appropriately
* run tests and examples
