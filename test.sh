#!/bin/bash
#../build/src/generate_training data/1069623.3.fna,data/1069624.3.fna,data/1069625.3.fna data/1069623.3.fna.gff3,data/1069624.3.fna.gff3,data/1069625.3.fna.gff3 5,1,2 -c 1 -r 250 -o test_rgi -a rgi_tsv
../build/src/generate_training data/1069623.3.fna,data/1069624.3.fna,data/1069625.3.fna data/1069623.3.fna.gff3,data/1069624.3.fna.gff3,data/1069625.3.fna.gff3 5,1,2 -c 1 -r 150 -o test_gff -e HSXt -a gff
#../build/src/generate_training data/1069623.3.fna,data/1069624.3.fna,data/1069625.3.fna data/1069623.3.fna.gff3,data/1069624.3.fna.gff3,data/1069625.3.fna.gff3 5,1,2 -c 1 -r 150 -o test_gff -e HSXt
