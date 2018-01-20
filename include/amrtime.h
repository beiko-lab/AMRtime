#ifndef AMRTIME_H 
#define AMRTIME_H

#include <string>
#include <seqan/arg_parse.h>

// ===========================================================================
// Classes
// ===========================================================================

class MainOptions {
    // class to contain formatted input options
    public:
        std::string mode;
        
};

// ===========================================================================
// Functions
// ===========================================================================

seqan::ArgumentParser::ParseResult parse_main_args(MainOptions& options, 
                                                   int argc,
                                                   char** argv);

int generate_training(int argc, char *argv[]);
#endif // #ifndef AMRTIME_H
