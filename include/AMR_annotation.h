#ifndef AMR_ANNOTATION_H 
#define AMR_ANNOTATION_H 

#include<iostream>
#include<string>
class AMR_annotation {
    // class to hold ARO annotation from RGI GFF output
    public:
        std::string contig;
        std::string aro;
        uint32_t start;
        uint32_t end;
        char strand;
        bool operator== (const AMR_annotation &other);
        std::ostream operator<< (std::ostream &os);
};

#endif // #ifndef AMR_ANNOTATION_H
