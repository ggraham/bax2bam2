cmake_minimum_required(VERSION 3.16)
project(bax2bam)

include(FetchContent)
include(ExternalProject)

FetchContent_Declare(
  pbbam
  GIT_REPOSITORY https://github.com/PacificBiosciences/pbbam.git
  GIT_TAG v1.6.0
)

FetchContent_Declare(
  pbcopper
  GIT_REPOSITORY https://github.com/PacificBiosciences/pbcopper.git
  GIT_TAG v1.9.0
)

FetchContent_Declare(
  blasr_libcpp
  GIT_REPOSITORY https://github.com/PacificBiosciences/blasr_libcpp.git
  GIT_TAG 5.3.4
)

FetchContent_Declare(
  htslib
  URL https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2
  URL_HASH MD5=c55c73099e2c5d71b084c267a9f20258
)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
include_directories(${HDF5_INCLUDE_DIRS})
#include_directories('/usr/include/hdf5/serial')

# ExternalProject_Add(
#   htslib
#   URL https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2
#   URL_HASH MD5=c55c73099e2c5d71b084c267a9f20258
#   DOWNLOAD_DIR _deps/htslib-dl
#   SOURCE_DIR   _deps/htslib-src
#   BINARY_DIR   _deps/htslib-bin
#   INSTALL_DIR  _deps/htslib-ins
# )

FetchContent_MakeAvailable(blasr_libcpp)
FetchContent_MakeAvailable(pbcopper)
FetchContent_MakeAvailable(pbbam)
FetchContent_MakeAvailable(htslib)

message(STATUS "INCLUDE LOCATION" ${HDF5_INCLUDE_DIRS})
message(STATUS "version" ${HDF5_VERSION})
message(STATUS "DEFINITIONS" ${HDF5_DEFINITIONS})
message(STATUS "LIBRARIES" ${HDF5_LIBRARIES})
message(STATUS "HL_LIBRARIES" ${HDF5_HL_LIBRARIES})

add_executable(${PROJECT_NAME} src/main.cpp)
target_link_libraries(${PROJECT_NAME} ${HDF5_HL_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES} blasr_libcpp pbcopper pbbam htslib)