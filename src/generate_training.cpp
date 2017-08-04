#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan; 

int count_nucleotides(std::string fasta_fp); 

class AMR_annotation {
    public:
        std::string aro;
        uint32_t start;
        uint32_t end;
        char strand;
};

std::map<std::string, std::vector<AMR_annotation>> read_amr_annotations(std::string gff_fp);



int main(int argc, char *argv[]){

    int genome_size = count_nucleotides(argv[1]);
    std::cout << genome_size << std::endl;     
    
    read_amr_annotations(argv[2]);

}




std::map<std::string, std::vector<AMR_annotation>> read_amr_annotations(std::string gff_fp) {
    // parse the RGI generated GFF file into a
    // dictionary of format {contig: [Annotation, Annotation, ...]}
    // for quick access

    GffFileIn gffFileIn;
    if (!open(gffFileIn, toCString(gff_fp))) {
        std::cerr << "ERROR: Could not open file: " << gff_fp << std::endl;
        //return -1;
    }

    GffRecord gffRecord;
    std::string aro;
    
    // contig: [start: 34, stop: 42,  
    std::<std::string, std::vector<AMR_annotation>> annotation_lookup;
    //
    std::vector<AMR_annotation> annotations;

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
            
            AMR_annotation annot {aro, 
                                  gffRecord.beginPos, 
                                  gffRecord.endPos, 
                                  gffRecord.strand};

            if( annotation_lookup.count(gffRecord.ref) > 0){
                // append to vector
                annotation_lookup.emplace(gffRecord.ref,
                                          annotation_lookup[gffRecord.ref].push_back(annot))
            }
            else {
                // create new 
                annotation_lookup.emplace(gffRecord.ref, 
                                          std::vector<AMR_annotation>{annot});
            }
            
        }
        catch (ParseError e) {
            break;
        }
    }
   return annotation_lookup;
}


int count_nucleotides(std::string fasta_fp) {
    
    int nt_count = 0;

    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(fasta_fp))) {
        std::cerr << "ERROR: Could not open file: " << fasta_fp << std::endl;
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
