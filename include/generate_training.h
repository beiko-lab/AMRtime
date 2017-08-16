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
        
};

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

std::vector<AMR_annotation> read_amr_annotations(std::vector<std::string> gff_list);

std::vector<std::string> split(std::string str, char delimiter);

std::string prepare_metagenome(std::vector<std::string> genome_list,
                               std::vector<uint32_t> abundance_list,
                               std::string output_name);

uint32_t count_nucleotides(std::string fasta_fp); 

void create_labels(std::vector<AMR_annotation> annotations, std::string sam_fp,
                   std::string output_name);

uint32_t estimate_read_depth(std::string combined_genome_fp,
                             uint32_t coverage_fold,
                             uint32_t read_length);

#endif // #ifndef GENERATE_TRAINING_H_
