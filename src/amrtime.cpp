#include <iostream>
#include "AMRtimeConfig.h"

int main(int argc, char const ** argv){
    
    std::cout << "Version info:" << std::endl;
    std::cout << AMRtime_VERSION_MAJOR 
        << AMRtime_VERSION_MINOR 
        << AMRtime_VERSION_PATCH 
        << std::endl;



    return 0;
}
