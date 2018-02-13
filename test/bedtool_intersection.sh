#!/bin/bash
sam2bed < output_errFree.sam > output_errFree.bed
bedtools intersect -a output_errFree.bed -b amr.bed -wo -f 0.333 
awk -F$'\t' 'BEGIN {OFS = FS}; {print substr($(NF-6), 1, length($(NF-6))-2), $(NF-3),$NF}' intersection.tsv > clean_hits.tsv
