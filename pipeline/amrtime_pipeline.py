#!/usr/bin/env python
from sklearn import preprocessing
from sklearn import model_selection

import pandas as pd

import subprocess
import argparse
import os
import numpy as np

from amrtime import utils
from amrtime import parsers
from amrtime import encoding


# how do I make params that only evaluate once again?
RANDOM_STATE = 42

def generate_training_data(card):

    if not os.path.exists('training_data'):
        os.mkdir('training_data')

    if not os.path.exists('training_data/card_proteins.faa'):
        card.write_proteins('training_data/card_proteins.faa')

    if not os.path.exists('training_data/card_nucleotides.fna'):
        card.write_nucleoties('training_data/card_nucleotides.fna')

    # generate training data
    if not os.path.exists('training_data/metagenome.fq'):
        subprocess.check_call('./art_illumina -q -na -ef -sam -ss MSv3 -i training_data/card_nucleotides.fna -f 2 -l 250 -rs 42 -o training_data/metagenome', shell=True)

    # build labels
    if not os.path.exists('training_data/metagenome_labels.tsv'):
        subprocess.check_call("grep '^@gb' training_data/metagenome.fq > training_data/read_names", shell=True)

        labels_fh = open('training_data/metagenome_labels.tsv', 'w')
        with open('training_data/read_names') as fh:
            for line in fh:
                line = line.strip().replace('@', '')
                read_name = line
                contig_name = '-'.join(line.split('-')[:-1])

                split_line = line.split('|')
                aro = split_line[2]
                amr_name = '-'.join(split_line[3].split('-')[:-1])

                # as they are all straight from the canonical sequence
                # and can't partially overlap
                amr_cutoff = 'Perfect'
                amr_overlap = '250'

                labels_fh.write('\t'.join([read_name, contig_name, aro, amr_name,
                                           amr_cutoff, amr_overlap]) + '\n')
        labels_fh.close()

    return 'training_data/metagenome.fq', 'training_data/metagenome_labels.tsv'


if __name__ == '__main__':

    parser = utils.get_parser()
    args = parser.parse_args()

    # ugh having an issue where the same sequence is appearing multiple times
    # in the card dna sequences thus leading to a mismatch in sequence number
    # when I key on accession during encoding/
    # why are the same sequences in card.dna_sequences multiple times?
    # have I fucked up my code or is that a thing in the CARD database
    if args.mode == 'train':
        card = parsers.CARD(args.card_fp)
        dataset, labels = generate_training_data(card)

        # run diamond filter to get
        homology_encoding = encoding.Homology(dataset, 'training_data/card_proteins.faa')
        X_sim = homology_encoding.encode(card)
        X_dissim = homology_encoding.encode(card, dissimilarity=True)

        #tnf = encoding.TNF('training_data/metagenome.fq')
        #X_tnf = tnf.encode()

        #np.save('training_data/X', X)
        #X = np.load('training_data/X.npy')
        aros = parsers.prepare_labels(labels, card)
        np.save('training_data/y', aros)
        y = np.load('training_data/y.npy')
        print(y)

        le = preprocessing.LabelEncoder()
        le.fit(y)

        #X_sim = pd.read_pickle('test_diamond_norm').as_matrix()
        #X_dissim = pd.read_pickle('test_diamond_norm_dissim').as_matrix()

        print(y.shape)
        print(X_sim.shape)
        print(X_dissim.shape)
        #model_selection.cross_val_score()
        #cv = model_selection.StratifiedKFold(n_splits=5, shuffle=True)
        #for train_index, test_index in cv.split(X_sim, y):
        #    clf = GaussianNB()
        #    clf.fit(X_sim[train_index], y[train_index])
        #    score = clf.score(X_sim[test_index], y[test_index])
        #    print(score)

        #for train_index, test_index in cv.split(X_dissim, y):
        #    clf.fit(X_dissim[train_index], y[train_index])
        #    score = clf.score(X_dissim[test_index], y[test_index])
        #    print(score)



        # rebalance training set
        # 5-fold CV (not bother with test-train split as we have other test-set)

        # create a dict to store gene family model classes keyed by the gene family
        #family_models = {}
        #for ix, aro in enumerate(aros):
        #    gene_family = card.aro_to_gene_family[aro][0]
        #    if gene_family not in family_models:
        #        family_models.update({gene_family: {'X': [X[ix,:]],
        #                                            'y': [aro]}
        #                             })
        #    else:
        #        family_models[gene_family]['X'].append(X[ix,:])
        #        family_models[gene_family]['y'].append(aro)

        #print(family_models)

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



