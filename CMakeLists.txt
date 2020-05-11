cmake_minimum_required(VERSION 2.6)

### use a compiler supporting constexpr math
set(CMAKE_C_COMPILER "/usr/local/lib/gcc-9.1.0/bin/gcc" CACHE STRING "gcc compiler" FORCE)
set(CMAKE_CXX_COMPILER "/usr/local/lib/gcc-9.1.0/bin/g++" CACHE STRING "g++ compiler" FORCE)

### Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

project(ShapedSlp)

include_directories(${PROJECT_SOURCE_DIR})
set(FOLCA_SOURCE_DIR ${PROJECT_SOURCE_DIR}/folca)
include_directories(${FOLCA_SOURCE_DIR})
set(SUX_SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/sux/sux)
include_directories(${SUX_SOURCE_DIR}/function)
include_directories(${SUX_SOURCE_DIR}/support)

### check for SDSL
find_library(SDSL_LIB sdsl)
if (NOT SDSL_LIB)
  message("sdsl-lite not found, building from submodule (hopefully you did git clone --recursive)")
  add_subdirectory(external/sdsl-lite)
  include_directories(${CMAKE_BINARY_DIR}/external/sdsl-lite/include/)
  include_directories(${CMAKE_BINARY_DIR}/external/sdsl-lite/external/libdivsufsort/include)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/external/sdsl-lite/lib)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/external/sdsl-lite/external/libdivsufsort/lib)
else()
  find_library(DSS divsufsort)
  find_library(DSS64 divsufsort64)
  message("sdsl: " ${SDSL_LIB})
  message("divsufsort: " ${DSS})
  message("divsufsort64: " ${DSS64})
endif()

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS "--std=c++2a /usr/lib64/libstdc++.so.6.0.27")
set(CMAKE_CXX_FLAGS "--std=c++2a")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -march=native -pthread")
set(CMAKE_CXX_FLAGS_RELEASE "-g -Ofast -fstrict-aliasing -march=native -pthread -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native -pthread")


add_executable(SlpEncBuild SlpEncBuild.cpp Common.cpp ${FOLCA_SOURCE_DIR}/fbtree.cpp ${FOLCA_SOURCE_DIR}/onlinebp.c)
TARGET_LINK_LIBRARIES(SlpEncBuild sdsl)

add_executable(SubstrBenchmark SubstrBenchmark.cpp Common.cpp ${FOLCA_SOURCE_DIR}/fbtree.cpp ${FOLCA_SOURCE_DIR}/onlinebp.c)
TARGET_LINK_LIBRARIES(SubstrBenchmark sdsl)
