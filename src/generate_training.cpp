#include <fstream>
#include <iostream>
#include <string>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan; 

int count_nucleotides(std::string fasta_fp); 

int main(int argc, char *argv[]){

    int genome_size = count_nucleotides(argv[1]);
    std::cout << genome_size << std::endl;     

    

}

int count_nucleotides(std::string fasta_fp){
    
    int nt_count = 0;

    SeqFileIn seqFileIn;
    if (!open(seqFileIn, toCString(fasta_fp)))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        return -1;
    }

    CharString id;
    Dna5String seq;
    try
    {
        while (!atEnd(seqFileIn))
        {
            readRecord(id, seq, seqFileIn);
            nt_count += length(seq);
        }
    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return -1;
    }

    return nt_count;
}
