#include <iostream>
#include "AMRtimeConfig.h"
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>


int main(int argc, char const ** argv){
    
    std::cout << "Version info:" << std::endl;

    std::cout << AMRtime_VERSION_MAJOR << "."
              << AMRtime_VERSION_MINOR << "." 
              << AMRtime_VERSION_PATCH << "."
              << std::endl;

    std::cout << seqan::CharString("Test") << std::endl;

    return 0;
}
