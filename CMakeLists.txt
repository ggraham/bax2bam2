cmake_minimum_required(VERSION 3.16)
project(bax2bam)

include(FetchContent)

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

FetchContent_MakeAvailable(blasr_libcpp)
FetchContent_MakeAvailable(pbcopper)
FetchContent_MakeAvailable(pbbam)