#!/usr/bin/env python

import subprocess
import argparse
import os

import utils
import parsers

# how do I make params that only evaluate once again?
RANDOM_STATE = 42

def generate_training_data(card_proteins)

    if not os.path.exists('training_data'):
        os.mkdir('training_data')

    synthetic_metagenome_fp = 'training_data/card_proteins.fasta'
    with open(synthetic_metagenome_fp, 'w') as fh:
        for accession, sequence in card_proteins:
            fh.write('>{}\n{}\n'.format(accession, sequence))

    # generate training data
    subprocess.check_output('art_illumina -q -na -ef -sam -ss MSv3 -i {} -f 10 -l 250 -rs 42 -o training_data/metagenome'.format(synthetic_metagenome_fp), shell=True)

    # build labels
    subprocess.check_output("grep '^@gb' training_data/metagenome.fq > training_data/read_names")

    labels_fh = open('training_data/metagenome_labels.tsv', 'w')
    with open('training_data/read_names') as fh:
        for line in fh:
            line = line.strip().replace('@', '')
            read_name = line
            contig_name = '-'.join(line.split('-')[:-1])

            split_line = line.split('|')
            aro = split_line[4].replace('ARO:', '')
            amr_name = '-'.join(split_line[5].split('-')[:-1])

            # as they are all straight from the canonical sequence
            # and can't partially overlap
            amr_cutoff = 'Perfect'
            amr_overlap = '250'

            labels_fh.write('\t'.join([read_name, contig_name, aro, amr_name,
                                       amr_cutoff, amr_overlap]) + '\n')
    labels_fh.close()
    print('Labels generated')


def prepare_dataset(dataset_fp):
    from parsers import MultiKModel

    filepath = 'dna2vec-20161219-0153-k3to8-100d-10c-29320Mbp-sliding-Xat.w2v'
    mk_model = MultiKModel(filepath)

    with open(dataset_fp) as fh:
        for ix, line in enumerate(fh):
            if ix % 4 == 1:
                seq = Read(line.strip(), mk_model)





if __name__ == '__main__':

    parser = utils.get_parser()
    args = parser.parse_args()

    if args.mode == 'train':
        card = parsers.CARD(args.card_fp)
        #dataset, labels = generate_data(card.proteins)
        dataset, labels = 'training_data/metagenome.fq', 'training_data/metagenome_labels.tsv'
        encode(dataset)
        # encode reads (feature extraction/kmer)
        # rebalance training set
        # 5-fold CV (not bother with test-train split as we have other test-set)
        # create a dict to store gene family model classes keyed by the gene family
        # for each gene family in map:
            # parse labels for indices of reads with aros in that family
            # create new x-y for those indices and encoded reads (out-of-core
            #   data-structure? dask?)
            # train model with 5-fold CV
            # add to dictionary keyed by gene family
    elif args.mode == 'test':
        pass
        # if read passes DIAMOND filter
        # encode read to kmers/features (maybe do in bulk for efficiency later)
        # run against trained gene family predictor
        # run against gene family model for highest predictor
        # what is max(y_pred) for read
        # evaluate against that gene family level model in the trained model
        #   dictionary



