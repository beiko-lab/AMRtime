# configures building lib

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
add_executable(generate_training generate_training.cpp)
target_link_libraries(generate_training ${SEQAN_LIBRARIES})
