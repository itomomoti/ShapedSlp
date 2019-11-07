cmake_minimum_required(VERSION 2.6)

### set a compiler supporting constexpr math
set(CMAKE_C_COMPILER "/opt/local/bin/gcc-mp-9" CACHE STRING "gcc compiler" FORCE)
set(CMAKE_CXX_COMPILER "/opt/local/bin/gcc-mp-9" CACHE STRING "g++ compiler" FORCE)

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

set(CMAKE_CXX_FLAGS "--std=c++2a /opt/local/lib/libstdc++.6.dylib")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-g -Ofast -fstrict-aliasing -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")


# add_executable(ShapedSlp ShapedSlp.cpp Common.cpp ${SUX_SOURCE_DIR}/function/DoubleEF.cpp ${SUX_SOURCE_DIR}/function/RiceBitVector.cpp ${SUX_SOURCE_DIR}/support/SpookyV2.cpp)
# TARGET_LINK_LIBRARIES(ShapedSlp sdsl)

# add_executable(ShapedSlp_SeparatedVlc ShapedSlp_SeparatedVlc.cpp Common.cpp ${SUX_SOURCE_DIR}/function/DoubleEF.cpp ${SUX_SOURCE_DIR}/function/RiceBitVector.cpp ${SUX_SOURCE_DIR}/support/SpookyV2.cpp)
# TARGET_LINK_LIBRARIES(ShapedSlp_SeparatedVlc sdsl)


# add_executable(ShapedSlp_Build ShapedSlp_Build.cpp Common.cpp ${SUX_SOURCE_DIR}/function/DoubleEF.cpp ${SUX_SOURCE_DIR}/function/RiceBitVector.cpp ${SUX_SOURCE_DIR}/support/SpookyV2.cpp)
# TARGET_LINK_LIBRARIES(ShapedSlp_Build sdsl)

# add_executable(ShapedSlp_Substr_Benchmark ShapedSlp_Substr_Benchmark.cpp Common.cpp ${SUX_SOURCE_DIR}/function/DoubleEF.cpp ${SUX_SOURCE_DIR}/function/RiceBitVector.cpp ${SUX_SOURCE_DIR}/support/SpookyV2.cpp)
# TARGET_LINK_LIBRARIES(ShapedSlp_Substr_Benchmark sdsl)


# add_executable(PoSlp_Build PoSlp_Build.cpp Common.cpp ${FOLCA_SOURCE_DIR}/fbtree.cpp ${FOLCA_SOURCE_DIR}/onlinebp.c)
# TARGET_LINK_LIBRARIES(PoSlp_Build sdsl)

# add_executable(PoSlp_Substr_Benchmark PoSlp_Substr_Benchmark.cpp Common.cpp ${FOLCA_SOURCE_DIR}/fbtree.cpp ${FOLCA_SOURCE_DIR}/onlinebp.c)
# TARGET_LINK_LIBRARIES(PoSlp_Substr_Benchmark sdsl)


add_executable(SlpEncBuild SlpEncBuild.cpp Common.cpp ${FOLCA_SOURCE_DIR}/fbtree.cpp ${FOLCA_SOURCE_DIR}/onlinebp.c ${SUX_SOURCE_DIR}/function/DoubleEF.cpp ${SUX_SOURCE_DIR}/function/RiceBitVector.cpp ${SUX_SOURCE_DIR}/support/SpookyV2.cpp)
TARGET_LINK_LIBRARIES(SlpEncBuild sdsl)

add_executable(SubstrBenchmark SubstrBenchmark.cpp Common.cpp ${FOLCA_SOURCE_DIR}/fbtree.cpp ${FOLCA_SOURCE_DIR}/onlinebp.c ${SUX_SOURCE_DIR}/function/DoubleEF.cpp ${SUX_SOURCE_DIR}/function/RiceBitVector.cpp ${SUX_SOURCE_DIR}/support/SpookyV2.cpp)
TARGET_LINK_LIBRARIES(SubstrBenchmark sdsl)
