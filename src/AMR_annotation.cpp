#include <iostream>
#include "AMR_annotation.h"

// custom comparison for AMR_Annotation that check for equality of all members
bool AMR_Annotation::operator== (const AMR_annotation &other ) const {
        std::vector<bool> comparison = {contig == other.contig,
                                           aro == other.aro,
                                         start == other.start,
                                           end == other.end,
                                        strand == other.strand};
        return comparison.all();
};

// Overload the output stream operator for annotation class
std::ostream AMR_Annotation::operator<< (std::ostream &os){
    // to dump the attributes
    return os << " contig: " << contig << " aro: " << aro 
            << " start: " << start << " end: " << end 
            << " strand: " << strand << std::endl;
};
