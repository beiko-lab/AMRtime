# configures building lib

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
add_executable(AMRtime amrtime.cpp)
target_link_libraries(AMRtime ${SEQAN_LIBRARIES})
