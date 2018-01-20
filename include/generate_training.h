#ifndef GENERATE_TRAINING_H_
#define GENERATE_TRAINING_H_

#include <seqan/arg_parse.h>
#include <vector>
#include <iostream>
#include <string>
#include <map>

// ===========================================================================
// Classes and Typedefs
// ===========================================================================
typedef std::vector<std::string> TStrList;

class GenerateOptions {
    // class to contain formatted input options
    public:
        TStrList genomes;
        std::vector<std::uint32_t> relative_abundances;
        TStrList annotations;
        uint32_t coverage;
        uint32_t read_length;
        std::string output_name;
        std::string annotation_type;
        std::string art_error_profile;
        uint32_t minimum_overlap;
};

// currently not used but will need later possible
class RgiRecord {
    // class to contain a record from the RGI TSV format
    //ORF_ID	Contig	Start	Stop	Orientation	Cut_Off	Pass_Bitscore	Best_Hit_Bitscore	Best_Hit_ARO	Best_Identities	ARO	ARO_name	Model_type	SNP	Best_Hit_ARO_category	ARO_category	Other_hit_bitscores	Predicted_DNA	Predicted_Protein	CARD_Protein_Sequence	ID	Model_ID

    public:
        std::string orf_id;
        std::string contig;
        uint32_t start;
        uint32_t stop;
        char orientation;
        std::string cut_off;
        double pass_bitscore;
        double best_hit_bitscore;
        std::string best_hit_aro;
        double best_identities;
        TStrList aro;
        TStrList aro_name;
        std::string model_type;
        std::string snp_best_hit_aro_category;
        std::string aro_category; 
        TStrList other_hit_bitscores;
        std::string predicted_DNA;
        std::string predicted_protein;
        std::string card_protein_sequence;
        std::string id;
        uint32_t model_id;
};


class AmrAnnotation {
    // class to hold ARO annotation from RGI GFF output
    public:
        std::string contig;
        std::string aro;
        std::string amr_name;
        std::string cutoff;
        uint32_t start;
        uint32_t end;
        char strand;
};

typedef std::map<std::string, std::vector<AmrAnnotation>> TAnnotationMap;

// ===========================================================================
// Functions
// ===========================================================================

seqan::ArgumentParser::ParseResult parseGenerateArgs(GenerateOptions& options, 
                                                     int argc,
                                                     char** argv);

int32_t rangeOverlap(uint32_t annot_start, uint32_t annot_end, 
                     uint32_t read_loc_start, uint32_t read_loc_end);

TAnnotationMap readAmrAnnotations(TStrList annotation_fps,
                                  std::string annotation_type);

TStrList split(std::string str, char delimiter);

std::string prepareMetagenome(TStrList genome_list,
                               std::vector<uint32_t> abundance_list,
                               std::string output_name);

void createLabels(TAnnotationMap annotations, 
                   std::string sam_fp,
                   std::string output_name,
                   uint32_t minimum_overlap);

uint32_t stoui32(const std::string& s);

int generateTraining(int argc, char *argv[]);


std::ostream& operator<< (std::ostream &out, const AmrAnnotation &annotation);
bool operator== (AmrAnnotation &first, const AmrAnnotation &other );
#endif // #ifndef GENERATE_TRAINING_H_
