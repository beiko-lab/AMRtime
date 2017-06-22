#include <iostream>
#include <sstream>
#include <string>
#include <seqan/arg_parse.h>
#include <unistd.h>
#include "AMRtimeConfig.h"

int main(int argc, char const ** argv){

    seqan::ArgumentParser parser("AMRtime");

    setShortDescription(parser, "Rapid AMR prediction from metagenomic datasets");

    setVersion(parser, AMRtime_VERSION);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIIN\fP \\fIDB\fP");

    addDescription(
            parser,
            "AMRtime is a tool to predict antimicrobial resistance (AMR) genes"
            "present in metagenomic 2nd-generation sequencing datasets in a "
            "rapid and reliable way.");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "IN"));

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "DB"));
  

    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK){
        return res == seqan::ArgumentParser::PARSE_ERROR;
    }
    
    seqan::CharString inputFile;
    getArgumentValue(inputFile, parser, 0);

    seqan::CharString dbFile;
    getArgumentValue(dbFile, parser, 1);

    std::stringstream ss;
    ss << "./diamond blastx -d " << dbFile << " -q " << inputFile << " -f 6 --min-score 60";
    std::cout << ss.str();
    system( ss.str().c_str() );
     
    
    return 0;
}
