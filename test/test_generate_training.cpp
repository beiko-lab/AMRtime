#include <gtest/gtest.h>
#include <stdexcept>
#include <iostream>
#include <cstdio>
#include "generate_training.h"

// ===========================================================================
// Utility Functions
// ===========================================================================

bool compareFiles(const std::string& p1, const std::string& p2) {
  std::ifstream f1(p1, std::ifstream::binary|std::ifstream::ate);
  std::ifstream f2(p2, std::ifstream::binary|std::ifstream::ate);

  if (f1.fail() || f2.fail()) {
    return false; //file problem
  }

  if (f1.tellg() != f2.tellg()) {
    return false; //size mismatch
  }

  //seek back to beginning and use std::equal to compare contents
  f1.seekg(0, std::ifstream::beg);
  f2.seekg(0, std::ifstream::beg);
  return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.rdbuf()));
}

// ===========================================================================
// Tests
// ===========================================================================

// Tests split correctly breaks string
TEST(splitTest, expected) {
    TStrList test_vec = {"A","B","C"};
    EXPECT_EQ(test_vec, split("A,B,C", ','));
    TStrList test_vec2 = {"123","Bas"};
    EXPECT_EQ(test_vec2, split("123;Bas", ';'));
    TStrList test_vec3 = {"#@","-2"};
    EXPECT_EQ(test_vec3, split("#@\t-2", '\t'));
}

// Tests range_overlap calculates the correct overlap size between ranges
TEST(rangeOverlapTest, expected){
    EXPECT_EQ(4, rangeOverlap(0, 19, 15, 39));
    EXPECT_EQ(0, rangeOverlap(0, 99, 250, 259));
    EXPECT_EQ(99, rangeOverlap(0, 499, 250, 349));
}

// Test string to uint32_t works correctly
TEST(stoui32Test, expected){
    EXPECT_EQ(15, stoui32("15"));
    EXPECT_THROW(stoui32("-1"), std::invalid_argument);
    EXPECT_EQ(4294967295,stoui32("4294967295")); 
    // Expect but don't really want
    EXPECT_EQ(4294967295, stoui32("9999999999"));
}

// Test readAmrAnnotation reads rgi tsv's correctly
TEST(readAmrAnnotation, rgi_tsv){

    TStrList test_annotation_tsvs = {"test/data/test_1.tsv", 
                                     "test/data/test_2.tsv"};

    // fill expected AMR annotation map
    TAnnotationMap expected;
    expected["AILI01000002"] = std::vector<AmrAnnotation> {
        AmrAnnotation {"AILI01000002", 
                       "3000616",
                       "mel", 
                       "Strict", 
                       60044, 
                       61507, 
                       '+'},
        AmrAnnotation {"AILI01000002", 
                       "3000822", 
                       "pmrA", 
                       "Perfect", 
                       74688, 
                       75887, 
                       '+'}
    };
    expected["ALCK01000005"] = std::vector<AmrAnnotation> {
        AmrAnnotation {"ALCK01000005", 
                       "3000822", 
                       "pmrA", 
                       "Strict", 
                       53713, 
                       54912, 
                       '+'},
        AmrAnnotation {"ALCK01000005", 
                       "3000025", 
                       "patB",
                       "Strict",
                       516527,
                       517693,
                       '-'}
    };

    TAnnotationMap actual = readAmrAnnotations(test_annotation_tsvs, "rgi_tsv");
    EXPECT_EQ(expected["AILI01000002"], actual["AILI01000002"]);
    EXPECT_EQ(expected["ALCK01000005"], actual["ALCK01000005"]);
}


// Test that prepareMetagenome makes the right file with the right coverage
TEST(prepareMetagenome, Correct){
    TStrList genome_list = {"test/data/test_1.fna", 
                            "test/data/test_2.fna"};
    
    std::string expected_output = "test/data/expected_metagenome.fna";
    std::string actual_output = prepareMetagenome(genome_list, 
                                                  std::vector<uint32_t> {3, 1},
                                                  "actual");
    EXPECT_TRUE(compareFiles(actual_output, expected_output));
    remove(actual_output.c_str());
}


// Test label generation (although only in the positive i.e. all the labels
// that should exist, existing)
TEST(createLabels, Correct){
    std::string sam_fp = "test/data/test.sam";
    TAnnotationMap annotations;

    annotations["AILI01000002"] = std::vector<AmrAnnotation> {
        AmrAnnotation {"AILI01000002", 
                       "3000616",
                       "mel", 
                       "Strict", 
                       60044, 
                       61507, 
                       '+'},

        AmrAnnotation {"AILI01000002", 
                       "3000822", 
                       "pmrA", 
                       "Perfect", 
                       74688, 
                       75887, 
                       '+'}
    };
    annotations["ALCK01000005"] = std::vector<AmrAnnotation> {
        AmrAnnotation {"ALCK01000005", 
                       "3000822", 
                       "pmrA", 
                       "Strict", 
                       53713, 
                       54912, 
                       '+'},
        AmrAnnotation {"ALCK01000005", 
                       "3000025", 
                       "patB",
                       "Strict",
                       516527,
                       517693,
                       '-'}
    };   
    
    std::string output = "test";
    createLabels(annotations, sam_fp, output, 50);
        
    EXPECT_TRUE(compareFiles("test_labels.tsv", 
                "test/data/bedtool_labels.tsv"));
    remove("test_labels.tsv");
}

// Test label generation (although only in the positive i.e. all the labels
// that should exist, existing)
TEST(getCleanReads, Correct){
    std::string actual_clean_fq = "test/data/bedtool_clean.fq";
    std::string actual_clean_labels = "test/data/bedtool_clean_labels.tsv";
    getCleanReads("test/data/bedtool");

    EXPECT_TRUE(compareFiles(actual_clean_fq,
                            "test/data/expected_clean.fq"));
    EXPECT_TRUE(compareFiles(actual_clean_labels,
                             "test/data/expected_clean_labels.tsv"));

    remove(actual_clean_labels.c_str());
    remove(actual_clean_fq.c_str());
}

// ./bin/amrtime generate_training test/data/test_1.fna,test/data/test_2.fna test/data/test_1.tsv,test/data/test_2.tsv 3,1 -x

// ===========================================================================
// Test Runner
// ===========================================================================

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
 
