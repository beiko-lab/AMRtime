#include "AMRtimeConfig.h"
#include "amrtime.h"
#include "generate_training.h"

#include <seqan/arg_parse.h>

// ===========================================================================
// Functions
// ===========================================================================

// ---------------------------------------------------------------------------
// Function parseMainArgs()
// ---------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult parseMainArgs(MainOptions& options, 
                                                 int argc,
                                                 char** argv){
    seqan::ArgumentParser parser("AMRtime");

    setShortDescription(parser, "Rapid AMR prediction from metagenomic datasets");

    setVersion(parser, AMRtime_VERSION);

    addUsageLine(parser, "\\fIMODE\\fP");

    addDescription(
            parser,
            "AMRtime is a tool to predict antimicrobial resistance (AMR) genes"
            "present in metagenomic 2nd-generation sequencing datasets in a "
            "rapid and reliable way.");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "mode"));
    setValidValues(parser, 0, "generate_training");

    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;
    }

    getArgumentValue(options.mode, parser, 0);
    return seqan::ArgumentParser::PARSE_OK;
}

// ==========================================================================
// Main 
// ==========================================================================

int main(int argc, char *argv[]){

    // TODO fix it so that the main argparse doesn't hijack internal ones
    /*seqan::ArgumentParser::ParseResult res = parseMainArgs(options, 
                                                           argc, 
                                                           argv);

    if (res != seqan::ArgumentParser::PARSE_OK){
        return res == seqan::ArgumentParser::PARSE_ERROR;
    }
    */

    // cut_off first arg and pass to mode 
    if(argc > 1){
        std::string mode = argv[1];
        if(mode == "generate_training"){
            int ret = generateTraining(argc - 1, argv + 1);
            return ret;
        } else {
            std::cout << "Mode invalid: must be one of [generate_training]" 
                << std::endl;
            return -1;
        }
    } else {
        std::cout << "Mode not specified: must be one of [generate_training]" 
            << std::endl;
        return -1;
    }
}
