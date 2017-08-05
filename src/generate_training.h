#ifndef GENERATE_TRAINING_H_
#define GENERATE_TRAINING_H_

class AMR_annotation {
    // class to hold ARO annotation from RGI GFF output
    public:
        std::string contig;
        std::string aro;
        uint32_t start;
        uint32_t end;
        char strand;
};


std::vector<AMR_annotation> read_amr_annotations(std::vector<std::string> gff_list);

std::vector<std::string> split(std::string str, char delimiter);

std::string prepare_metagenome(std::vector<std::string> genome_list,
                               std::vector<uint32_t> abundance_list,
                               std::string output_name);

int count_nucleotides(std::string fasta_fp); 

void create_labels(std::vector<AMR_annotation> annotations, std::string sam_fp);

uint32_t estimate_read_depth(std::string combined_genome_fp,
                             uint32_t coverage_fold,
                             uint32_t read_length);

#endif // #ifndef GENERATE_TRAINING_H_
