#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include "generate_training.h"

using namespace seqan; 

int main(int argc, char *argv[]){
    // generate_training genome1.fasta,genome2.fasta,genome3.fasta 4,2,3 genome1.gff,genome2.gff,genome3.gff

    std::vector<std::string> genome_list = split(argv[1], ',');
    std::vector<std::string> abundance_list = split(argv[2], ',');

    std::vector<uint32_t> abundances;
    for (uint32_t i=0; i<length(abundance_list); ++i) {
        int abundance = atoi(abundance_list.at(i).c_str());
        abundances.push_back(abundance);
    }

    prepare_metagenome(genome_list, abundances);
    
    std::vector<std::string> gff_list = split(argv[3], ',');
    std::vector<AMR_annotation> amr_annotations = read_amr_annotations(gff_list);
    return 0;
}

std::vector<AMR_annotation> read_amr_annotations(std::vector<std::string> gff_list) {

    std::vector<AMR_annotation> annotations;

    for(uint32_t i = 0; i < length(gff_list); ++i){
        GffFileIn gffFileIn;
        if (!open(gffFileIn, toCString(gff_list[i]))) {
            std::cerr << "ERROR: Could not open file: " << gff_list[i] << std::endl;
            //return -1;
        }

        GffRecord gffRecord;
        std::string aro;

        while (!atEnd(gffFileIn)) {
            try {
                readRecord(gffRecord, gffFileIn);

                aro = "";
                
                for (unsigned i = 0; i < length(gffRecord.tagValues[1]); ++i) {
                    if (gffRecord.tagValues[1][i] == ','){
                        break;
                    }
                    aro.push_back(gffRecord.tagValues[1][i]);
                }
                
                AMR_annotation annotation {toCString(gffRecord.ref),
                                           aro, 
                                           gffRecord.beginPos, 
                                           gffRecord.endPos, 
                                           gffRecord.strand};

                annotations.push_back(annotation);

            }
            catch (ParseError e) {
                break;
            }
        }
   }
   return annotations;
}

int prepare_metagenome(std::vector<std::string> genome_list,
                       std::vector<uint32_t> abundance_list) {
    // Copy the genomes up to necessary numbers into the artifical
    // metagenome contigs
    
    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;
    
    // read each genome in list and 'amplify' for relative abundance
    for(uint32_t i = 0; i < length(genome_list); ++i){
        
        SeqFileIn seqFileIn;
        if (!open(seqFileIn, toCString(genome_list[i]))) {
            std::cerr << "ERROR: Could not open file: " << genome_list[i] << std::endl;
            return -1;
        }

        StringSet<CharString> temp_ids;
        StringSet<Dna5String> temp_seqs;

        readRecords(temp_ids, temp_seqs, seqFileIn);
        
        // append the sequences as many times as the relative abundance
        // implies i.e. if its 3, copy the sequences 3x into the master
        // metagenome fasta
        for(uint32_t j = 0; j < abundance_list[i]; ++j){
            append(ids, temp_ids);
            append(seqs, temp_seqs);
        }

    }
    
    // dump artificial metagenome to single fasta file
    SeqFileOut seqFileOut;
    if (!open(seqFileOut, "temp_metagenome.fasta")) {
            std::cerr << "ERROR: Could not open file: temp_metagenome.fasta" << std::endl;
            return -1;
    }

    writeRecords(seqFileOut, ids, seqs);

    return 0;
    
}
    
int count_nucleotides(std::string combined_genome_fp){
    // read a fasta and count the nucleotides
    int nt_count = 0;
    
    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(combined_genome_fp))) {
        std::cerr << "ERROR: Could not open file: " << combined_genome_fp << std::endl;
        return -1;
    }

    CharString id;
    Dna5String seq;
    while (!atEnd(seqFileIn)){
        readRecord(id, seq, seqFileIn);
        nt_count += length(seq);
    }

    return nt_count;
}


std::vector<std::string> split(std::string str, char delimiter) {
    // split a string on a specific delimiter
    std::vector<std::string> split_string;
    std::stringstream ss(str);
    std::string fragment;

    while(std::getline(ss, fragment, delimiter)){
        split_string.push_back(fragment);
    }
    
    return split_string;
}
