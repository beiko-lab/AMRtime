# configures dependencies, platform specifics and output paths

cmake_minimum_required(VERSION 3.7.2)
project(AMRtime)

# versioning
set(AMRtime_VERSION_MAJOR 0)
set(AMRtime_VERSION_MINOR 0)
set(AMRtime_VERSION_PATCH 1)
set(AMRtime_VERSION "${AMRtime_VERSION_MAJOR}.${AMRtime_VERSION_MINOR}.${AMRtime_VERSION_PATCH}")

# include_dirs
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set(PROJECT_TEST_DIR "${PROJECT_SOURCE_DIR}/test")
include_directories(${PROJECT_INCLUDE_DIR})

# config header
configure_file (
    "${PROJECT_INCLUDE_DIR}/AMRtimeConfig.h.in"
    "${PROJECT_BINARY_DIR}/AMRtimeConfig.h"
    )
include_directories("${PROJECT_BINARY_DIR}")


# dependencies
find_package(ZLIB)

find_package(SeqAn REQUIRED)
include_directories(${SEQAN_INCLUDE_DIRS})
add_definitions(${SEQAN_DEFINITIONS})

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# add subdirectory cmakelists
add_subdirectory(src)
#add_subdirectory(test)
