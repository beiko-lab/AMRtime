#-*- coding: utf-8 -*-
#!/usr/bin/env python

import argparse
import subprocess
import os

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prototype of AMRtime \
            metagenomic antimicrobial resistance analysis tool')

    parser.add_argument('input', type=str, help='Path to input fastq file')
    parser.add_argument('db', type=str, help='Path to CARD protein database')
    parser.add_argument('minscore', type=int, help='Minimum bitscore')

    args = parser.parse_args()

    output = args.input + '.out'
    cmd = "./diamond blastx --db {} -q {} -f 6 -p 2 \
            --min-score {} -o {} ".format(args.db,
                                          args.input,
                                          args.minscore,
                                          output)

    subprocess.call(cmd, shell=True)

    temp = output + '.temp'
    if os.path.exists(temp):
        os.remove(temp)

    print('cutting')
    cmd = 'cut -f 1 {} > {}'.format(output, temp)
    subprocess.call(cmd, shell=True)


    filtered_input = args.input + '.{}_filtered'.format(args.minscore)

    print('filtering')
    cmd = './seqtk subseq {} {} > {}'.format(args.input,
                                             temp,
                                             filtered_input)
    subprocess.call(cmd, shell=True)
