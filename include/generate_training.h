#ifndef GENERATE_TRAINING_H_
#define GENERATE_TRAINING_H_

class Options {
    // class to contain formatted input options
    public:
        std::vector<std::string> genomes;
        std::vector<std::uint32_t> relative_abundances;
        std::vector<std::string> annotations;
        uint32_t coverage;
        uint32_t read_length;
        std::string output_name;
        std::string annotation_type;
        std::string art_error_profile;
        
};

class RGI_record {
    // class to contain a record from the RGI TSV format
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
        double best_identity;
        std::vector<std::string> aro;
        std::vector<std::string> aro_name;
        std::string model_type;
        std::string snp;
        std::vector<std::string> best_hit_aro_category;
        std::string aro_category; 
        std::vector<std::string> other_hit_bitscores;
        std::string predicted_DNA;
        std::string predicted_protein;
        std::string card_protein_sequence;
        std::string id;
        uint32_t model_id;
};


     //ORF_ID	Contig	Start	Stop	Orientation	Cut_Off	Pass_Bitscore	Best_Hit_Bitscore	Best_Hit_ARO	Best_Identities	ARO	ARO_name	Model_type	SNP	Best_Hit_ARO_category	ARO_category	Other_hit_bitscores	Predicted_DNA	Predicted_Protein	CARD_Protein_Sequence	ID	Model_ID

class AMR_annotation {
    // class to hold ARO annotation from RGI GFF output
    public:
        std::string contig;
        std::string aro;
        uint32_t start;
        uint32_t end;
        char strand;
};

int32_t range_overlap(uint32_t annot_start, uint32_t annot_end, 
                       uint32_t read_loc_start, uint32_t read_loc_end);

std::vector<AMR_annotation> read_amr_annotations(std::vector<std::string> gff_list,
                                                 std::string annotation_type);

std::vector<std::string> split(std::string str, char delimiter);

std::string prepare_metagenome(std::vector<std::string> genome_list,
                               std::vector<uint32_t> abundance_list,
                               std::string output_name);

uint32_t count_nucleotides(std::string fasta_fp); 

void create_labels(std::vector<AMR_annotation> annotations, std::string sam_fp,
                   std::string output_name);


#endif // #ifndef GENERATE_TRAINING_H_
