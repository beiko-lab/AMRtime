# configures building lib
add_executable(amrtime amrtime.cpp generate_training.cpp)
target_link_libraries(amrtime ${SEQAN_LIBRARIES})
