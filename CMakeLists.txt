# 
# CMakeLists.txt : DilateAF
# Copyright (C) 2018  Prakhar Jaiswal
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

