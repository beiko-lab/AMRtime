#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

#include "generate_training.h"

using namespace seqan;

struct Options {
    std::vector<std::string> genomes;
    std::vector<std::uint32_t> relative_abundances;
    std::vector<std::string> annotations;
    uint32_t coverage;
    uint32_t read_length;
    std::string output_name;
}

seqan::ArgumentParser::ParseResult parse_command_line(Options & options, 
                                                      int argc,
                                                      char const ** argv){
    
// generate_training genome1.fasta,genome2.fasta,genome3.fasta 4,2,3 genome1.gff,genome2.gff,genome3.gff coverageX read_length output_name
    //
    //
    seqan::ArgumentParser parser("generate_training");

    setShortDescription(parser, "Synthetic Metagenomes Generator");

    setVersion(parser, AMRtime_VERSION);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIGENOME_LIST\fP \\fIANNOTATION_LIST\fP \\fIABUNDANCE_LIST\fP");

    addDescription(
            parser,
            "Tool to generate synthetic metagenomes at specified coverage"
            " and relative abundances from annotated genomes.");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "genomes"));

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "annotations"));

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "abundances"));

    addOption(parser, seqan::ArgParseOption(
        seqan::ArgParseArgument::INT, "coverage",
        "Required coverage for metagenome"));
    setDefaultValue(parser, "coverage", 1);

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INT, "read_length",
        "length of reads to simulate"));
    setDefaultValue(parser, "read_length", 150);

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "output_name",
        "output file name"));
    setDefaultValue(parser, "output_name", "output");

    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;
    }
    
    std::string temp;
    getArgumentValue(temp, parser, "genomes");
    options.genomes = split(temp);

    getArgumentValue(temp, parser, "annotations");
    options.annotations = split(temp);

    getArgumentValue(temp, parser, "abundances");
    std::vector<std::string> abundance_strings = split(temp);
    for (uint32_t i=0; i<abundance_strings.size(); ++i) {
        uint32_t abundance = std::stoi(abundance_strings.at(i).c_str());
        options.relative_abundances.push_back(abundance);
    } 

    bool length_ok = (options.genomes.size() == options.annotations.size() && \
            options.annotations.size () == options.relative_abundances.size());
    if !(length_ok) {
        std::cerr << "ERROR: You must provide the same number of genomes, annotations (and relative abundances if specified)\n";
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(options.coverage, parser, "coverage");
    getOptionValue(options.read_length, parser, "read_length");
    getOptionValue(options.output_name, parser, "output_name");

    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char *argv[]){

    Options options;
    seqan::ArgumentParser::ParseResult res = parse_command_line(options, 
                                                                argc, 
                                                                argv)

    if (res != seqan::ArgumentParser::PARSE_OK){
        return res == seqan::ArgumentParser::PARSE_ERROR;
    }
    
    
    std::vector<std::string> genome_list = split(argv[1], ',');
    std::vector<std::string> abundance_list = split(argv[2], ',');

    std::uint32_t coverage_fold = std::stoi(argv[3]);
    std::uint32_t read_length = std::stoi(argv[4]);
    std::string output_name = argv[5];

    std::vector<uint32_t> abundances;
    for (uint32_t i=0; i<length(abundance_list); ++i) {
        int abundance = std::stoi(abundance_list.at(i).c_str());
        abundances.push_back(abundance);
    }


    std::string metagenome_fp = prepare_metagenome(genome_list, abundances,
                                                   output_name);


    uint32_t read_number = estimate_read_depth(metagenome_fp,
                                     coverage_fold,
                                     read_length);

    std::stringstream ss;
    std::string simulated_sam_fp = output_name + ".sam";
    ss << "mason_simulator -ir " << metagenome_fp << " -n " << read_number << " -oa " << simulated_sam_fp << " -o " << output_name + ".fq" << std::endl;
    std::cout << ss.str() << std::endl;
    system(ss.str().c_str());

    std::vector<std::string> gff_list = split(argv[3], ',');
    std::vector<AMR_annotation> amr_annotations = read_amr_annotations(gff_list);

    create_labels(amr_annotations, simulated_sam_fp);

    return 0;
}


void create_labels(std::vector<AMR_annotation> annotations, std::string sam_fp){

    BamFileIn samFileIn(sam_fp.c_str());
    BamAlignmentRecord record;
    while (!atEnd(samFileIn)){
        readRecord(record, samFileIn);
        // check annotations
        // if in any range of simulated read print list of AROs to label file
        // else print 0 on line
    }
}

std::vector<AMR_annotation> read_amr_annotations(std::vector<std::string> gff_list) {

    std::vector<AMR_annotation> annotations;

    for(uint32_t i = 0; i < length(gff_list); ++i){
        GffFileIn gffFileIn;
        if (!open(gffFileIn, toCString(gff_list[i]))) {
            std::cerr << "ERROR: Could not open file: " << gff_list[i] << std::endl;
            exit(1);
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

std::string prepare_metagenome(std::vector<std::string> genome_list,
                               std::vector<uint32_t> abundance_list,
                               std::string output_name){
    // Copy the genomes up to necessary numbers into the artifical
    // metagenome contigs

    StringSet<CharString> ids;
    StringSet<Dna5String> seqs;

    // read each genome in list and 'amplify' for relative abundance
    for(uint32_t i = 0; i < length(genome_list); ++i){

        SeqFileIn seqFileIn;
        if (!open(seqFileIn, toCString(genome_list[i]))) {
            std::cerr << "ERROR: Could not open file: " << genome_list[i] << std::endl;
            exit(1);
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
    std::string metagenome_fp = output_name + "metagenome.fasta";
    SeqFileOut seqFileOut;
    if (!open(seqFileOut, toCString(metagenome_fp))) {
            std::cerr << "ERROR: Could not open file: temp_metagenome.fasta" << std::endl;
            exit(1);
    }

    writeRecords(seqFileOut, ids, seqs);

    return metagenome_fp;

}

uint32_t count_nucleotides(std::string combined_genome_fp){
    // read a fasta and count the nucleotides
    int nt_count = 0;

    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(combined_genome_fp))) {
        std::cerr << "ERROR: Could not open file: " << combined_genome_fp << std::endl;
        exit(1);
    }

    CharString id;
    Dna5String seq;
    while (!atEnd(seqFileIn)){
        readRecord(id, seq, seqFileIn);
        nt_count += length(seq);
    }

    return nt_count;
}


uint32_t estimate_read_depth(std::string combined_genome_fp,
                             uint32_t coverage_fold,
                             uint32_t read_length){
    // Calculate the necessary read number for a specified read depth

    uint32_t nt_count = count_nucleotides(combined_genome_fp);


    double read_number_float = (coverage_fold * nt_count) / read_length;

    uint32_t read_number = static_cast<int>(read_number_float);

    return read_number;

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
