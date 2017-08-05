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

int prepare_metagenome(std::vector<std::string> genome_list, std::vector<uint32_t> abundance_list);

int count_nucleotides(std::string fasta_fp); 

#endif // #ifndef GENERATE_TRAINING_H_
