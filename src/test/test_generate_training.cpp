#include<gtest/gtest.h>
#include <stdexcept>
#include "generate_training.h"

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
    EXPECT_EQ(5, rangeOverlap(0, 19, 15, 39));
    EXPECT_EQ(0, rangeOverlap(0, 99, 250, 259));
    EXPECT_EQ(100, rangeOverlap(0, 499, 250, 349));
}

// Test string to uint32_t works correctly
TEST(stoui32Test, expected){
    EXPECT_EQ(15, stoui32("15"));
    EXPECT_THROW(stoui32("-1"), std::invalid_argument);
    EXPECT_EQ(4294967295,stoui32("4294967295")); 
    // Expect but don't really want
    EXPECT_EQ(4294967295, stoui32("9999999999"));
}

//TAnnotationMap readAmrAnnotations(TStrList annotation_fps,
//                                 std::string annotation_type);
TEST(readAmrAnnotation, rgi_tsv){
    
    TStrList test_annotation_tsvs = {"test_data/test_1.tsv", 
                                     "test_data/test_2.tsv"};
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



//std::string prepareMetagenome(TStrList genome_list,
//                               std::vector<uint32_t> abundance_list,
//                               std::string output_name);
TEST(prepareMetagenome, Correct){
    EXPECT_EQ(1,1);
}

//void createLabels(TAnnotationMap annotations, 
//                  std::string sam_fp,
//                  std::string output_name);
TEST(createLabels, Correct){
    EXPECT_EQ(1,1);
}



int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
 
