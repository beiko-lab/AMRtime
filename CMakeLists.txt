cmake_minimum_required (VERSION 3.7.2)
project (AMRtime)

# versioning
set (AMRtime_VERSION_MAJOR 0)
set (AMRtime_VERSION_MINOR 0)
set (AMRtime_VERSION_PATCH 1)

configure_file (
    "${PROJECT_SOURCE_DIR}/include/AMRtimeConfig.h.in"
    "${PROJECT_BINARY_DIR}/AMRtimeConfig.h"
    )

include_directories("${PROJECT_BINARY_DIR}")

# add subdirectory cmakelists
add_subdirectory(src)
add_subdirectory(test)
