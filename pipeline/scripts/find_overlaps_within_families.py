#!/usr/bin/env python

import glob
import itertools
import re
from tqdm import tqdm
import pickle

def longest_common_substring(s1, s2):
   m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
   longest, x_longest = 0, 0
   for x in range(1, 1 + len(s1)):
       for y in range(1, 1 + len(s2)):
           if s1[x - 1] == s2[y - 1]:
               m[x][y] = m[x - 1][y - 1] + 1
               if m[x][y] > longest:
                   longest = m[x][y]
                   x_longest = x
           else:
               m[x][y] = 0
   return s1[x_longest - longest: x_longest]


if __name__ == '__main__':

    lcs = {}
    #for family in glob.glob('analyses/family_fasta/*'):
    family = 'analyses/TEM_beta-lactamase_clean_unique.fna'
    family_seqs = []
    with open(family) as fh:
        ix = 0
        for line in fh:
            if ix % 2 == 0:
                acc = line.strip()
                name = acc.split('|')[3]
            elif ix % 2 == 1:
                seq = line.strip()
                family_seqs.append((acc, seq))
            ix+=1

    lcs_lens = []
    for subset in tqdm(itertools.combinations(family_seqs, 2)):
        lcs = longest_common_substring(subset[0][1], subset[1][1])
        lcs_lens.append((subset[0][0], subset[1][0], len(lcs)))

    with open('tem_lens', 'wb') as fh:
        pickle.dump(lcs_lens, fh)

