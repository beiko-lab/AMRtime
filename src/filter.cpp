#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <stdexcept>

#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/arg_parse.h>
#include <seqan/sequence.h>

#include "AMRtimeConfig.h"
#include "generate_training.h"

#define RANDOM_SEED 42

// ===========================================================================
// Functions
// ===========================================================================

// ---------------------------------------------------------------------------
// Function parseGenerateArgs()
// ---------------------------------------------------------------------------
seqan::ArgumentParser::ParseResult parseGenerateArgs(GenerateOptions& options, 
                                                     int argc,
                                                     char** argv){
    /* Build argument parser, validate and store parsed arguments 
     *in an instance of the GenerateOptions class.
     */
    seqan::ArgumentParser parser("filter");

    setShortDescription(parser, "Metagenomic Filtering");

    setVersion(parser, AMRtime_VERSION);

    setDate(parser, "August 2018");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fImteagnome.fq\\fP \\fICARD.json\\fP");

    addDescription(
            parser,
            "Tool to filter a metagenome for possible AMR determinants before"
            " sensitive classification.");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "data_tsv"));

    addOption(parser, seqan::ArgParseOption(
        "a", "annotation_type", "File type of RGI annotations.",
        seqan::ArgParseOption::STRING, "annotation_type"));
    setValidValues(parser, "annotation_type", "gff rgi_tsv");
    setDefaultValue(parser, "annotation_type", "rgi_tsv");

    addOption(parser, seqan::ArgParseOption(
        "e", "art_error_profile", "Read error profile.",
        seqan::ArgParseOption::STRING, "art_error_profile"));
    setValidValues(parser, "art_error_profile", "MSv3 HS25");
    setDefaultValue(parser, "art_error_profile", "HS25");

    addOption(parser, seqan::ArgParseOption(
        "c", "coverage", "Required coverage for metagenome.",
        seqan::ArgParseOption::INTEGER, "coverage"));
    setMinValue(parser, "coverage", "1");
    setDefaultValue(parser, "coverage", 1);
    
    // as 150 is the max length for HiSeq error profile again
    // should be a dependency in the limits for theses
    addOption(parser, seqan::ArgParseOption(
        "r", "read_length", "Length of reads to simulate.",
        seqan::ArgParseOption::INTEGER, "read_length"));
    setMinValue(parser, "read_length", "33");
    setDefaultValue(parser, "read_length", 150);

    addOption(parser, seqan::ArgParseOption(
        "o", "output_name", "Output file name.",
        seqan::ArgParseOption::STRING, "output_name"));
    setDefaultValue(parser, "output_name", "output");
    
    // in theory read length should be the maximum which is currently not set
    addOption(parser, seqan::ArgParseOption(
        "m", "minimum_overlap", "Minimum annotation overlap to label.",
        seqan::ArgParseOption::INTEGER, "min_overlap"));
    setMinValue(parser, "minimum_overlap", "1");
    setDefaultValue(parser, "minimum_overlap", 50);

    addOption(parser, seqan::ArgParseOption(
        "x", "get_clean_reads", 
        "Optionally output only positive labelled data"));

    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;
    }
    
    
    // parse data specified in tsv
    std::string data_fp;
    getArgumentValue(data_fp, parser, 0);
    TStrList split_line;
    std::string line;
    std::ifstream data_fh (data_fp);
    
    // read data file 
    if (data_fh.is_open()){
        while (std::getline(data_fh, line)){
                split_line = split(line, '\t');

                // check for valid 3 column format
                if(split_line.size() != 3){
                    std::cerr << "Data file must be 3 column tsv formatted " <<
                        data_fp << std::endl;
                    std::exit(1);
                }
                options.genomes.push_back(split_line[0]);
                options.annotations.push_back(split_line[1]);
                options.relative_abundances.push_back(std::stoi(split_line[2].c_str()));

        }
    } else {
        std::cerr << "ERROR: Could not open metadata file: " << data_fp 
            << std::endl;
       std::exit(1);
    }

    bool length_ok = (options.genomes.size() == options.annotations.size() && \
            options.annotations.size () == options.relative_abundances.size());
    if (!length_ok) {
        std::cerr << "ERROR: You must provide the same number "
                     "of genomes, annotations and relative abundances "
                     "if specified)" << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(options.coverage, parser, "coverage");
    getOptionValue(options.read_length, parser, "read_length");
    getOptionValue(options.output_name, parser, "output_name");
    getOptionValue(options.art_error_profile, parser, "art_error_profile");
    getOptionValue(options.annotation_type, parser, "annotation_type");
    getOptionValue(options.minimum_overlap, parser, "minimum_overlap");
    getOptionValue(options.get_clean_reads, parser, "get_clean_reads");
    return seqan::ArgumentParser::PARSE_OK;
}


