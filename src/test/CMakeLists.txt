enable_testing()
add_executable(run_unit_tests test_generate_training.cpp ../generate_training.cpp)

configure_file(test.sh test.sh COPYONLY)
configure_file(test_data/expected_metagenome.fna 
               test_data/expected_metagenome.fna 
               COPYONLY)
configure_file(test_data/test_1.fna test_data/test_1.fna COPYONLY)
configure_file(test_data/test_1.tsv test_data/test_1.tsv COPYONLY)
configure_file(test_data/test_2.fna test_data/test_2.fna COPYONLY)
configure_file(test_data/test_2.tsv test_data/test_2.tsv COPYONLY)

target_link_libraries(run_unit_tests GTest::Main ${SEQAN_LIBRARIES})
add_test(NAME run_unit_tests COMMAND run_unit_tests)
