#
# cpp_api/CMakeLists.txt
#
#
# The MIT License
#
# Copyright (c) 2018 TileDB, Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

cmake_minimum_required(VERSION 3.3)
project(DilateAF
  VERSION 3.5.0
  LANGUAGES CXX)

find_package(ArrayFire)
find_package(TileDB)

# Set C++11 as required standard for all C++ targets.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (WIN32)
  # Disable MSVC warnings not present with gcc.
  #   C4204: non-constant aggregate initializer
  #   C4221: aggregate initializer with local variable
  #   C4305: implicit truncation of double to float
  add_compile_options(/wd4204 /wd4221 /wd4305)
else()
  # Don't treat printf format warnings as errors.
  add_compile_options(-Wno-error=format)
endif()

add_executable(dilateAF ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/foam.cpp)

if(ArrayFire_CUDA_FOUND)
    target_link_libraries(dilateAF ArrayFire::afcuda)
elseif(ArrayFire_OpenCL_FOUND)
    target_link_libraries(dilateAF ArrayFire::afopencl)
elseif(ArrayFire_CPU_FOUND)
    target_link_libraries(dilateAF ArrayFire::afcpu)
endif()

target_include_directories(dilateAF PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/)
target_link_libraries(dilateAF TileDB::tiledb_shared)

# On Linux, must explicitly link -lpthread -ldl in order for static linking
# to libzstd or libcurl to work.
target_link_libraries(dilateAF pthread dl)

