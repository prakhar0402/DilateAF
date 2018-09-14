# Copyright (c) 2018, ArrayFire
# All rights reserved.
#
# This file is distributed under 3-clause BSD license.
# The complete license agreement can be obtained at:
# http://arrayfire.com/licenses/BSD-3-Clause

cmake_minimum_required(VERSION 3.0)
project(DilateAF
  VERSION 3.5.0
  LANGUAGES CXX)

find_package(ArrayFire)

if(ArrayFire_CPU_FOUND)
  add_executable(dilateAF_cpu src/main.cpp)
  target_link_libraries(dilateAF_cpu ArrayFire::afcpu)
endif()

if(ArrayFire_CUDA_FOUND)
  add_executable(dilateAF_cuda src/main.cpp)
  target_link_libraries(dilateAF_cuda ArrayFire::afcuda)
endif()

if(ArrayFire_OpenCL_FOUND)
  add_executable(dilateAF_opencl src/main.cpp)
  target_link_libraries(dilateAF_opencl ArrayFire::afopencl)
endif()
