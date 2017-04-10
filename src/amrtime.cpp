#include <iostream>
#include <seqan/arg_parse.h>
#include "AMRtimeConfig.h"

int main(int argc, char const ** argv){

    seqan::ArgumentParser parser("AMRtime");

    setCategory(parser, "AMRtime");

    setShortDescription(parser, "Rapid AMR prediction from metagenomic datasets");

    setVersion(parser, AMRtime_VERSION);

    addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fIIN\fP");

    addDescription(
            parser,
            "AMRtime is a tool to predict antimicrobial resistance (AMR) genes"
            "present in metagenomic 2nd-generation sequencing datasets in a "
            "rapid and reliable way.");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "IN"));
  

    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK){
        return res == seqan::ArgumentParser::PARSE_ERROR;
    }
    
    seqan::CharString inputFile;
    getArgumentValue(inputFile, parser, 0);
    
    std::cout << inputFile << std::endl;

    return 0;
}
