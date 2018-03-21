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
    /* Custom ostream operator to dump all the attribute to output
     * for the AmrAnnotation class
     */
    return out << "contig: " << annotation.contig 
            << ", aro: " << annotation.aro 
            << ", amr_name: " << annotation.amr_name 
            << ", cutoff: " << annotation.cutoff
            << ", start: " << annotation.start 
            << ", end: " << annotation.end 
            << ", strand: " << annotation.strand;
};

bool operator== (const AmrAnnotation &first, const AmrAnnotation &other) {
    /* Custom comparison for the AmrAnnotation class that just checks
     * all the attributes match
     */
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
    } else {
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
    /* Build argument parser, validate and store parsed arguments 
     *in an instance of the GenerateOptions class.
     */
    seqan::ArgumentParser parser("generate_training");

    setShortDescription(parser, "Synthetic Metagenomes Generator");

    setVersion(parser, AMRtime_VERSION);

    setDate(parser, "January 2018");

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIData_File.tsv\\fP");

    addDescription(
            parser,
            "Tool to generate synthetic metagenomes at specified coverage"
            " and relative abundances from annotated genomes as specified"
            " in a supplied tsv formatted as genome_fp\tannotation_fp\tabundance.");

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
        std::getline(data_fh, line);
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

// ---------------------------------------------------------------------------
// Function split()
// ---------------------------------------------------------------------------
TStrList split(std::string str, char delimiter) {
    /* split a string on a specific delimiter
     */
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
    /* Copy the genomes up to necessary numbers into the artifical
     * metagenome contigs
     */

    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;

    // open filehandle to artificial metagenome to dump single fasta file
    std::string metagenome_fp = output_name + "_synthetic_metagenome.fasta";
    seqan::SeqFileOut seqFileOut (metagenome_fp.c_str());

    // read each genome in list and 'amplify' for relative abundance
    for(uint32_t genome_ix = 0; genome_ix < genome_list.size(); ++genome_ix){
        
        // if you use the constructor it gives appropriate error
        seqan::SeqFileIn seqFileIn(genome_list[genome_ix].c_str());

        seqan::StringSet<seqan::CharString> temp_ids;
        seqan::StringSet<seqan::Dna5String> temp_seqs;

        readRecords(temp_ids, temp_seqs, seqFileIn);

        // append the sequences as many times as the relative abundance
        // implies i.e. if its 3, copy the sequences 3x into the master
        // metagenome fasta
        // append the copy number to the id as a suffix to prevent sam parsing 
        // errors
        for(uint32_t copy_number = 0; 
                copy_number < abundance_list[genome_ix]; 
                ++copy_number){
            
            // copy the sequences in the specified number of times
            // add the correct copy number to accession
            // add it to the first space as SAM etc truncate at the first
            // space of the accession
            seqan::StringSet<seqan::CharString> ids_with_copy_number;
            for(auto &id: temp_ids){
                std::stringstream id_ss;
                id_ss << id;
                std::string id_str = id_ss.str();

                std::size_t first_space = id_str.find_first_of(" ");
                std::string first_part = id_str.substr(0, first_space);
                std::string second_part = id_str.substr(first_space);

                std::stringstream id_with_copy;
                id_with_copy << first_part << "_" << copy_number << second_part;

                appendValue(ids_with_copy_number, id_with_copy.str());
            }

            writeRecords(seqFileOut, ids_with_copy_number, temp_seqs);
        }

    }

    return metagenome_fp;

}

// ---------------------------------------------------------------------------
// Function readAmrAnnotations()
// ---------------------------------------------------------------------------
TAnnotationMap readAmrAnnotations(
        TStrList annotation_fps, 
        std::string annotation_type) {
    /* Parse the RGI TSV or RGI GFF file and store the AMR annotation details
     * in a contig-keyed map with lists of AmrAnnotation classes as the values.
     */

    std::vector<AmrAnnotation> annotations;
    
    // if gff files then use the gff parser
    if(annotation_type == "gff"){
        for (auto &gff_fp : annotation_fps){

            seqan::GffFileIn gffFileIn (gff_fp.c_str());
            seqan::GffRecord gff_record;
            std::string aro;

            while (!atEnd(gffFileIn)) {
                try {
                    readRecord(gff_record, gffFileIn);

                    aro = "";
                    
                    for (uint32_t i = 0; 
                            i < length(gff_record.tagValues[1]); 
                            ++i){
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
                    
                    //remove just the copy number from the contig name
                    std::string contig_name = split_line[1].substr(0, 
                                              split_line[1].find_last_of("_"));

                    AmrAnnotation annotation {contig_name,
                                              split_line[10],
                                              split_line[11],
                                              split_line[5],
                                              stoui32(split_line[2]),
                                              stoui32(split_line[3]),
                                              split_line[4].c_str()[0]};
                    annotations.push_back(annotation);
                }
            } else {
                std::cerr << "ERROR: Could not open file: " 
                    << tsv_fp << std::endl;
                std::exit(1);
            }
        }
   } else {
         std::cerr << "ERROR: invalid annotation type: " 
             << annotation_type << std::endl;
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
    /* Convert string to unsigned 32-bit integer
     * make sure it isn't negative for some reason as this shouldn't 
     * happen in the input as they are coords.
     */
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
    /* Using the parsed annotations, check the locations for the simulated
     * reads in the SAM file and determine whether they overlap with the 
     * annotated AMR genes
     * 
     * Output the labels in the format of a tsv 
     * read_name \t contig_name (and which copy) \t aro \t name \t cut-off \t overlap
     */
    
    // Open input file, BamFileIn can read SAM and BAM files.
    seqan::BamFileIn bamFileIn (sam_fp.c_str());

    // Open output file, BamFileOut accepts also an ostream and a format tag.
    std::ofstream labels_fh;
    labels_fh.open (output_name + "_labels.tsv");

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
            // only check when the annotation's contig is the same as the reads
            std::stringstream contig_name_ss;
            contig_name_ss << seqan::toCString(contigNames(bamContext)[bam_record.rID]);
            std::string contig_name = contig_name_ss.str();
            
            // strip out contig copy number indicator
            std::string contig_without_copy = contig_name.substr(0, 
                                                contig_name.find_last_of("_"));

            for (auto &annotation : annotations[contig_without_copy]){
                int32_t overlap = rangeOverlap(annotation.start, 
                                               annotation.end,
                                               bam_record.beginPos, 
                                               bam_record.beginPos + length(bam_record.seq));

                if (overlap > minimum_overlap) {
                    std::stringstream label_ss;
                        label_ss << bam_record.qName << "\t"
                            << contig_name << "\t"
                            << annotation.aro << "\t" 
                            << annotation.amr_name << "\t"
                            << annotation.cutoff << "\t" 
                            << overlap << std::endl;
                    labels.push_back(label_ss.str());
                }
            }

            // remove any potential duplicate labels, unecessary?
            std::sort(labels.begin(), labels.end());
            labels.erase(std::unique(labels.begin(), labels.end()), labels.end());
            
            // print the labels
            if (labels.size() == 0) {
                labels_fh << bam_record.qName << "\t"
                    << contig_name << "\t"
                    << "na" << "\t"
                    << "na" << "\t"
                    << "na" << "\t"
                    << "na" << std::endl;
            }
            else {
                for(auto &label: labels){
                    labels_fh << label;
                }
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
    /* Calculate the overlap size between two sets of start and end coordinates
     */

    int minimum = std::max(annot_start, read_loc_start);
    int maximum = std::min(annot_end, read_loc_end);

    if (minimum <= maximum) {
        // adding 1 to make it an inclusive range overlap
        return maximum - minimum;
    }
    else {
        //std::cout << 0 << std::endl;
        return 0;
    }
}

// ---------------------------------------------------------------------------
// Function getCleanReads()
// ---------------------------------------------------------------------------
void getCleanReads(std::string output_name){
    /* Simultaneously read through the label file and output fastq
     * and dump them into a new file iff there is a positive label (i.e.
     * an ARO for that read.
     */
    
    // open all the required input and output files
    std::string output_fq = output_name + ".fq";
    seqan::SeqFileIn seqFileIn(output_fq.c_str());

    std::string clean_fq = output_name + "_clean.fq";
    seqan::SeqFileOut seqFileOut(clean_fq.c_str());
    
    std::string output_label_fp = output_name + "_labels.tsv";
    std::ifstream output_labels(output_label_fp);
    
    std::ofstream clean_labels;
    clean_labels.open(output_name + "_clean_labels.tsv");

    // initialise variables to hold the values as I'm reading over the files
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::CharString qual;

    std::string label;
    std::string overlap;
    std::string id_string;
   
    // get first read
    readRecord(id, seq, qual, seqFileIn);
    id_string = seqan::toCString(id);

    // get first label
    getline(output_labels, label); 
    TStrList label_data = split(label, '\t');
       
    // loop over reads 
    while (!atEnd(seqFileIn)) {
        
        // loop over all annotations applying to the same read
        // i.e. the same read can have multiple annotations
        while(label_data[0] == id_string){

            // if actually labelled write to clean files
            if(label_data[2] != "na"){
                writeRecord(seqFileOut, id, seq, qual);
                clean_labels << label << std::endl;
            }
            
            
            // once we've done the first annotation grab the next one
            getline(output_labels, label); 
            label_data = split(label, '\t');
        }

        readRecord(id, seq, qual, seqFileIn);
        id_string = seqan::toCString(id);
    }
    clean_labels.close();
}


// ==========================================================================
// Function generateTraining()
// ==========================================================================
int generateTraining(int argc, char *argv[]){
    /* Main runner function for the generation of training data mode
     * parameterised by the post 'mode' selection command line args.
     *
     * It creates a synthetic metagenome and reads RGI output for the input
     * genomes to determine whether a given read in the synthetic metagenome
     * is from an AMR determinant.
     */

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

    std::cout << "Parsing Annotations: ";
    for(auto &annotation : options.annotations) {
        std::cout << annotation << " ";
    }
    std::cout << std::endl << std::endl;

    TAnnotationMap amr_annotations = readAmrAnnotations(options.annotations,
                                                        options.annotation_type);
    
    // get labels for error free for now
    std::cout << "Creating Labels: " 
        << options.output_name + "_labels.tsv " 
        << std::endl << std::endl;

    createLabels(amr_annotations, 
                 errfree_simulated_sam_fp, 
                 options.output_name, 
                 options.minimum_overlap);

    if(options.get_clean_reads){
        std::cout << "Creating Clean Labels: " 
            << options.output_name + "_clean.fq " 
            << options.output_name + "_clean_labels.tsv" << std::endl;
        getCleanReads(options.output_name);
    }

    return 0;
}
