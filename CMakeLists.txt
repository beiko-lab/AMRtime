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
set(SEQAN_FIND_DEPENDENCIES ZLIB BZip2)
find_package(SeqAn REQUIRED)

include_directories(${SEQAN_INCLUDE_DIRS})
add_definitions(${SEQAN_DEFINITIONS})

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# set CMAKE C++ flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -DSEQAN_DISABLE_VERSION_CHECK=YES")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# add subdirectory cmakelists
add_subdirectory(src)

# add tests subdirectoy
if (test)
  add_subdirectory(test)
endif()
