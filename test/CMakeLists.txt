enable_testing()
add_executable(run_unit_tests test_generate_training.cpp ../src/generate_training.cpp)

configure_file(test.sh test.sh COPYONLY)
configure_file(data/expected_metagenome.fna 
               data/expected_metagenome.fna 
               COPYONLY)
configure_file(data/test.sam data/test.sam COPYONLY)
configure_file(data/test_1.fna data/test_1.fna COPYONLY)
configure_file(data/test_1.tsv data/test_1.tsv COPYONLY)
configure_file(data/test_2.fna data/test_2.fna COPYONLY)
configure_file(data/test_2.tsv data/test_2.tsv COPYONLY)

configure_file(data/bedtool_labels.tsv data/bedtool_labels.tsv COPYONLY)
configure_file(data/bedtool.fq data/bedtool.fq COPYONLY)
configure_file(data/expected_clean.fq data/expected_clean.fq COPYONLY)
configure_file(data/expected_clean_labels.tsv data/expected_clean_labels.tsv COPYONLY)

target_link_libraries(run_unit_tests GTest::Main ${SEQAN_LIBRARIES})
add_test(NAME test COMMAND run_unit_tests)


