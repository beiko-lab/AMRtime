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
// Overwritten operators for class
// ===========================================================================

std::ostream& operator<< (std::ostream &out, const AmrAnnotation &annotation){
    // to dump the attributes
    return out << "contig: " << annotation.contig << ", aro: " << annotation.aro 
            << ", amr_name: " << annotation.amr_name << ", cutoff: " << annotation.cutoff
            << ", start: " << annotation.start << ", end: " << annotation.end 
            << ", strand: " << annotation.strand;
};

bool operator== (AmrAnnotation &first, const AmrAnnotation &other ) {
        bool comparison[7] = {first.contig == other.contig,
                              first.aro == other.aro,
                              first.amr_name == other.amr_name,
                              first.cutoff == other.cutoff,
                              first.start == other.start,
                              first.end == other.end,
                              first.strand == other.strand};
        // quick sum as comparison
        int sum = std::accumulate(comparison, comparison + 7, 0);

        // all true
        if(sum == 7) {
            return true;
        }
        else {
            return false;
        }
};

// ===========================================================================
// Functions
// ===========================================================================

// ---------------------------------------------------------------------------
// Function parseGenerateArgs()
// ---------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult parseGenerateArgs(GenerateOptions& options, 
                                                     int argc,
                                                     char** argv){
    // Build argument parser and store parsed arguments 
    // in an instance of the Options class
    seqan::ArgumentParser parser("generate_training");

    setShortDescription(parser, "Synthetic Metagenomes Generator");

    setVersion(parser, AMRtime_VERSION);

    setDate(parser, "January 2018");

    addUsageLine(parser, "[\\fIOPTIONS] \\fIGENOMELIST\f \\fIANNOTATIONLIST\f \\fIABUNDANCELIST\f");

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

    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;
    }
    
    std::string temp;
    getArgumentValue(temp, parser, 0);
    options.genomes = split(temp, ',');

    getArgumentValue(temp, parser, 1);
    options.annotations = split(temp, ',');

    getArgumentValue(temp, parser, 2);
    TStrList abundance_strings = split(temp, ',');
    for (uint32_t i=0; i<abundance_strings.size(); ++i) {
        uint32_t abundance = std::stoi(abundance_strings.at(i).c_str());
        options.relative_abundances.push_back(abundance);
    } 

    bool length_ok = (options.genomes.size() == options.annotations.size() && \
            options.annotations.size () == options.relative_abundances.size());
    if (!length_ok) {
        std::cerr << "ERROR: You must provide the same number of genomes, annotations "
                  << "(and relative abundances if specified)" << std::endl;
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    getOptionValue(options.coverage, parser, "coverage");
    getOptionValue(options.read_length, parser, "read_length");
    getOptionValue(options.output_name, parser, "output_name");
    getOptionValue(options.art_error_profile, parser, "art_error_profile");
    getOptionValue(options.annotation_type, parser, "annotation_type");
    getOptionValue(options.minimum_overlap, parser, "minimum_overlap");
    return seqan::ArgumentParser::PARSE_OK;
}

// ---------------------------------------------------------------------------
// Function split()
// ---------------------------------------------------------------------------

TStrList split(std::string str, char delimiter) {
    // split a string on a specific delimiter
    TStrList split_string;
    std::stringstream ss(str);
    std::string fragment;

    while(std::getline(ss, fragment, delimiter)){
        split_string.push_back(fragment);
    }

    return split_string;
}

// ---------------------------------------------------------------------------
// Function prepareMetagenome()
// ---------------------------------------------------------------------------

std::string prepareMetagenome(TStrList genome_list,
                               std::vector<uint32_t> abundance_list,
                               std::string output_name) {
    // Copy the genomes up to necessary numbers into the artifical
    // metagenome contigs

    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;

    // read each genome in list and 'amplify' for relative abundance
    for(uint32_t genome_ix = 0; genome_ix < genome_list.size(); ++genome_ix){

        seqan::SeqFileIn seqFileIn;
        if (!open(seqFileIn, seqan::toCString(genome_list[genome_ix]))) {
            std::cerr << "ERROR: Could not open file: " << genome_list[genome_ix] << std::endl;
            exit(1);
        }

        seqan::StringSet<seqan::CharString> temp_ids;
        seqan::StringSet<seqan::Dna5String> temp_seqs;

        readRecords(temp_ids, temp_seqs, seqFileIn);

        // append the sequences as many times as the relative abundance
        // implies i.e. if its 3, copy the sequences 3x into the master
        // metagenome fasta
        // append the copy number to the id as a suffix to prevent sam parsing errors
        for(uint32_t copy_number = 0; copy_number < abundance_list[genome_ix]; ++copy_number){
            append(ids, temp_ids);
            append(seqs, temp_seqs);
        }

    }

    // dump artificial metagenome to single fasta file
    std::string metagenome_fp = output_name + "metagenome.fasta";
    seqan::SeqFileOut seqFileOut;
    if (!open(seqFileOut, seqan::toCString(metagenome_fp))) {
            std::cerr << "ERROR: Could not open file: temp_metagenome.fasta" << std::endl;
            exit(1);
    }

    writeRecords(seqFileOut, ids, seqs);

    return metagenome_fp;

}

// ---------------------------------------------------------------------------
// Function readAmrAnnotations()
// ---------------------------------------------------------------------------
TAnnotationMap readAmrAnnotations(
        TStrList annotation_fps, 
        std::string annotation_type) {

    std::vector<AmrAnnotation> annotations;
    
    // if gff files then use the gff parser
    if(annotation_type == "gff"){
        for (auto &gff_fp : annotation_fps){

            seqan::GffFileIn gffFileIn;
            if (!open(gffFileIn, seqan::toCString(gff_fp))) {
                std::cerr << "ERROR: Could not open file: " << gff_fp << std::endl;
                std::exit(1);
            }

            seqan::GffRecord gff_record;
            std::string aro;

            while (!atEnd(gffFileIn)) {
                try {
                    readRecord(gff_record, gffFileIn);

                    aro = "";
                    
                    for (uint32_t i = 0; i < length(gff_record.tagValues[1]); ++i) {
                        if (gff_record.tagValues[1][i] == ','){
                            break;
                        }
                        aro.push_back(gff_record.tagValues[1][i]);
                    }
                    
                    // truncate the gff suffix from the contig name
                    std::string contig_name = seqan::toCString(gff_record.ref);
                    contig_name = contig_name.substr(0, contig_name.find("_"));
                    
                    // build annotation data together
                    // amr_name and cut-off not supported for this 
                    AmrAnnotation annotation {contig_name,
                                               aro,
                                               "unsupported", 
                                               "unsupported", 
                                               gff_record.beginPos,
                                               gff_record.endPos,
                                               gff_record.strand};

                    annotations.push_back(annotation);


                }
                // necessary for \r and \n endline chars?
                catch (seqan::ParseError const & e) {
                    break;
                }
            }
        }    
   }
   
   // if using the custom rgi TSV with the 
   else if(annotation_type == "rgi_tsv"){
        for (auto &tsv_fp : annotation_fps){

            std::string line;
            TStrList split_line;
            std::ifstream tsv_fh (tsv_fp);
            if (tsv_fh.is_open()){
                // skip header
                std::getline(tsv_fh, line);
                while (std::getline(tsv_fh, line)){

                split_line = split(line, '\t');
                    
                    AmrAnnotation annotation {split(split_line[1], '_')[0],
                                               split_line[10],
                                               split_line[8],
                                               split_line[5],
                                               stoui32(split_line[2]),
                                               stoui32(split_line[3]),
                                               split_line[4].c_str()[0]};
                    std::cout << annotation << std::endl;
                    annotations.push_back(annotation);
                }
            } else {
                std::cerr << "ERROR: Could not open file: " << tsv_fp << std::endl;
                std::exit(1);
            }
        }
   } else {
         std::cerr << "ERROR: invalid annotation type: " << annotation_type << std::endl;
          exit(1);
   }

   // make annotations way more efficient to query when creating labels
   // by first making a map for each contig to the amr set 
   TAnnotationMap annotation_map;
   for (auto &annotation : annotations){
       annotation_map[annotation.contig].push_back(annotation);
   }

   return annotation_map;
}

// ---------------------------------------------------------------------------
// Function stoui32()
// ---------------------------------------------------------------------------

uint32_t stoui32(const std::string& s) {
    // make sure it isn't negative for some reason
    if(s[0] == '-'){
        throw std::invalid_argument("Received negative value");
    }
    std::istringstream reader(s);
    uint32_t val = 0;
    reader >> val;
    return val;
}


// ---------------------------------------------------------------------------
// Function createLabels()
// ---------------------------------------------------------------------------

void createLabels(TAnnotationMap annotations, 
                  std::string sam_fp, 
                  std::string output_name,
                  uint32_t minimum_overlap){
    
    // Open input file, BamFileIn can read SAM and BAM files.
    seqan::BamFileIn bamFileIn;
    if (!open(bamFileIn, seqan::toCString(sam_fp))) {
        std::cerr << "ERROR: Could not open " << sam_fp << std::endl;
        std::exit(1);
    }

    // Open output file, BamFileOut accepts also an ostream and a format tag.
    std::ofstream labels_fh;
    labels_fh.open (output_name + ".labels");

    std::ofstream overlaps_fh;
    overlaps_fh.open (output_name + ".overlaps");
    
    try {
        // Copy header.
        seqan::BamHeader header;
        readHeader(header, bamFileIn);

        // read header into queriable context
        seqan::BamAlignmentRecord bam_record;
        typedef seqan::FormattedFileContext<seqan::BamFileIn, void>::Type TBamContext;
        TBamContext const & bamContext = context(bamFileIn);

        // for every read generated
        while (!atEnd(bamFileIn)) {
            readRecord(bam_record, bamFileIn);
            
            // create a vector of labels
            TStrList labels {};
            std::vector<uint32_t> overlaps {};
            
            // by checking all the annotations
            for (auto &annotation : annotations[seqan::toCString(contigNames(bamContext)[bam_record.rID])]){
                
                // only check when the annotation's contig is the same as the reads
                //if (annotation.contig == seqan::toCString(contigNames(bamContext)[bamRecord.rID])) {
                    bool both_pos = annotation.strand == '+' and not hasFlagRC(bam_record);
                    bool both_neg = annotation.strand == '-' and hasFlagRC(bam_record);
                    
                    // both on the same strand
                    if (both_pos or both_neg) {
                        
                        int32_t overlap = rangeOverlap(annotation.start, 
                                                       annotation.end,
                                                       bam_record.beginPos, 
                                                       bam_record.beginPos + length(bam_record.seq));
                        if (overlap > minimum_overlap) {
                            labels.push_back(annotation.aro);
                            overlaps.push_back(overlap); 
                        }
                    }
            }
            // sort and find unique AROs in labels due gff duplication issue wit RGI
            std::sort(labels.begin(), labels.end());
            labels.erase(std::unique(labels.begin(), labels.end()), labels.end());
            
            // print the labels
            if (labels.size() == 0) {
                labels_fh << "NONE" << std::endl;
            }
            if (labels.size() == 0) {
                overlaps_fh << "0" << std::endl;
            }
            else {
                for (auto &label : labels) {
                    labels_fh << label << ' ';
                }
                labels_fh << std::endl;
                
                for (auto &overlap : overlaps){
                    overlaps_fh << overlap << ' ';
                }
                overlaps_fh << std::endl;
            }
        }
    }
    
    
    catch (seqan::Exception const & e) {
        std::cout << "ERROR: " << e.what() << std::endl;
        std::exit(1);
    }
    
    labels_fh.close();
}

// ---------------------------------------------------------------------------
// Function rangeOverlap()
// ---------------------------------------------------------------------------

int32_t rangeOverlap(uint32_t annot_start, uint32_t annot_end, 
                     uint32_t read_loc_start, uint32_t read_loc_end){
    
    //std::cout << annot_start << " " << annot_end << std::endl;
    //std::cout << read_loc_start << " " << read_loc_end << std::endl;

    int minimum = std::max(annot_start, read_loc_start);
    int maximum = std::min(annot_end, read_loc_end);

    if (minimum <= maximum) {
        // adding 1 to make it an inclusive range overlap
        return maximum - minimum + 1;
    }
    else {
        //std::cout << 0 << std::endl;
        return 0;
    }
}

// ==========================================================================
// Function generateTraining()
// ==========================================================================

int generateTraining(int argc, char *argv[]){
    // Create simulated metagenome and label the reads

    // Get command line arguments
    GenerateOptions options;
    seqan::ArgumentParser::ParseResult res = parseGenerateArgs(options, 
                                                               argc, 
                                                               argv);

    if (res != seqan::ArgumentParser::PARSE_OK){
        return res == seqan::ArgumentParser::PARSE_ERROR;
    }
    
    // Generate the synthetic metagenome fasta file
    // i.e. concatenate the input genomes and 'amplify'
    // them up to the required copy number to satisfy relative
    // abundances
    std::cout << "Creating Synthetic Metagenome Fasta: ";
    for(auto &genome : options.genomes) {
        std::cout << genome << " ";
    }
    std::cout << std::endl << std::endl;

    std::string metagenome_fp = prepareMetagenome(options.genomes, 
                                                  options.relative_abundances,
                                                  options.output_name);

    
    // Simulate the reads themselves using mason and a system call
    std::cout << "Simulating Reads: " << options.read_length << "bp " 
        << options.coverage << "X coverage " << options.art_error_profile 
        << " error profile" << std::endl << std::endl;
    std::string simulated_sam_fp = options.output_name + ".sam";
    std::string errfree_simulated_sam_fp = options.output_name + "_errFree.sam";
    

    // run art simulator with and without MiSeq error profiles
    std::stringstream art_cmd;
    art_cmd << "art_illumina -q -na -ef -sam -ss " << options.art_error_profile 
        << " -i " << metagenome_fp << " -l " << options.read_length 
        << " -f " <<  options.coverage << " -o " << options.output_name 
        << " -rs " << RANDOM_SEED
        << " 2>&1 > /dev/null" << std::endl;
    system(art_cmd.str().c_str());

    //// get the error free reads using picard as art doesn't output
    // std::stringstream picard_cmd;
    //picard_cmd << "java -jar picard.jar SamToFastq I=" << errfree_simulated_sam_fp << 
    //    " FASTQ=" << options.output_name + "_errFree.fq 2> /dev/null" << std::endl;
    //System(picard_cmd.str().c_str());

    std::cout << "Parsing annotations: ";
    for(auto &annotation : options.annotations) {
        std::cout << annotation << " ";
    }
    std::cout << std::endl << std::endl;

    TAnnotationMap amr_annotations = readAmrAnnotations(options.annotations,
                                                        options.annotation_type);
    
    // get labels for error free for now
    std::cout << "Creating labels: " << options.output_name + ".labels" << std::endl;
    createLabels(amr_annotations, errfree_simulated_sam_fp, 
                 options.output_name, options.minimum_overlap);

    return 0;
}
