#!/usr/bin/env python
from sklearn import preprocessing
from sklearn import model_selection
from sklearn import ensemble, naive_bayes, tree, neighbors, neural_network, discriminant_analysis, linear_model, svm, gaussian_process
from sklearn import metrics

import pandas as pd

import subprocess
import argparse
import os
import numpy as np
from tqdm import tqdm

from amrtime import utils
from amrtime import parsers
from amrtime import encoding

# how do I make params that only evaluate once again?
RANDOM_STATE = 42

def classification_report_csv(report, fp):
    report_data = []
    lines = report.split('\n')
    for line in lines[2:-3]:
        row = {}
        row_data = line.split('      ')
        row['class'] = row_data[-5]
        row['precision'] = float(row_data[-4])
        row['recall'] = float(row_data[-3])
        row['f1_score'] = float(row_data[-2])
        row['support'] = float(row_data[-1])
        report_data.append(row)
    dataframe = pd.DataFrame.from_dict(report_data)
    #print(dataframe)
    dataframe.to_csv(fp, index = False)


def train_family_classifier(X_train, y_train):
    print("Training Family Classifier")
    clf = ensemble.RandomForestClassifier()
    clf.fit(X_train, y_train)
    return clf


def try_model(X, y, name):
    X_train, X_test, y_train, y_test = model_selection.train_test_split(
                                    X, y, stratify=y, test_size=0.33, random_state=42)

    clf = ensemble.RandomForestClassifier()
    try:
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)

        print(name, clf.__class__.__name__, metrics.precision_recall_fscore_support(y_test, y_pred, average='weighted'))
        classification_report_csv(metrics.classification_report(y_test, y_pred), name + clf.__class__.__name__)
    except:
        print(name, clf.__class__.__name__, ' failed')


def generate_training_data(card):

    if not os.path.exists('training_data'):
        os.mkdir('training_data')

    if not os.path.exists('training_data/card_proteins.faa'):
        card.write_proteins('training_data/card_proteins.faa')

    if not os.path.exists('training_data/card_nucleotides.fna'):
        card.write_nucleoties('training_data/card_nucleotides.fna')

    # generate training data
    if not os.path.exists('training_data/metagenome.fq'):
        subprocess.check_call('art_illumina -q -na -ef -sam -ss MSv3 -i training_data/card_nucleotides.fna -f 2 -l 250 -rs 42 -o training_data/metagenome', shell=True)

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


def train_family_level_classifiers(X_aro_train, y_aro_train, card):
    """
    Train the family level classifiers
    """
    print("Training ARO classifiers for each family")
    family_level_classifiers = {}
    for family in tqdm(X_family_bitscore_norm.columns):

        # get all the aros relevant to the family
        family_aros = card.gene_family_to_aro[family]

        # filter input to just the columns containing similarity to aros
        # within the family
        X_train = X_aro_train[family_aros]

        # get the indices where the label is one of the AROs belonging to
        # the current family
        label_indices = [ix for ix, x in enumerate(y_aro_train) if x in family_aros]

        y_train = np.array(y_aro_train)[label_indices]

        # grab only the reads where the label index is an ARO belonging to the
        # current family being trained
        X_train = X_train.iloc[label_indices]

        family_clf = ensemble.RandomForestClassifier()

        family_clf.fit(X_train, y_train)

        family_level_classifiers.update({family: [family_clf, family_aros]})

    return family_level_classifiers


def score(family_clf, family_classifiers, X_family, X_aro, y_family, y_aro):

        # predict and calculate performance for X
        print('Predicting Families')
        y_family_pred = family_clf.predict(X_family)
        classification_report_csv(metrics.classification_report(y_family,
                                                                y_family_pred),
                                  "family")

        y_aro_pred = []
        print('Gather predictions')
        # this takes way too long so instead of running individual preds
        # let's first group into families
        family_preds = {}
        for ix, family in enumerate(y_family_pred):
            if family in family_preds:
                family_preds[family].append(ix)
            else:
                family_preds.update({family: [ix]})

        print("Predicting ARO")
        all_indices = []
        for family in tqdm(family_preds):

            indices_for_preds_of_family = family_preds[family]

            all_indices += indices_for_preds_of_family
            family_level_clf = family_classifiers[family][0]
            family_aros = family_classifiers[family][1]

            X_family_sub = X_aro[family_aros]
            X_family_sub = X_family_sub.iloc[indices_for_preds_of_family]

            pred = family_level_clf.predict(X_family_sub)
            y_aro_pred.append(pred)

        y_aro_pred = np.hstack(y_aro_pred)

        classification_report_csv(metrics.classification_report(np.array(y_aro)[all_indices],
                                                                y_aro_pred),
                                  "aros")


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
        #homology_encoding = encoding.Homology(dataset, 'training_data/card_proteins.faa', 'DIAMOND')
        #X_family_bitscore_norm, X_aro_bitscore_norm = homology_encoding.encode(card, 'bitscore', norm=True)


        X_family_bitscore_norm = pd.read_pickle('family_bitscore.pkl')
        X_aro_bitscore_norm = pd.read_pickle('aro_bitscore.pkl')

        #X_pident = homology_encoding.encode(card, 'pident')
        #X_bitscore = homology_encoding.encode(card, 'bitscore')
        #X_evalue = homology_encoding.encode(card, 'evalue')
        #X_dissim = homology_encoding.encode(card, 'bitscore', dissimilarity=True)
        #X_dissim_norm = homology_encoding.encode(card, 'bitscore', norm=True, dissimilarity=True)
        #tnf = encoding.Kmer('training_data/metagenome.fq', 4)
        #X_tnf = tnf.encode()
        #k5mer = encoding.Kmer('training_data/metagenome.fq', 5)
        #X_5mer = k5mer.encode()

        amr_family_labels, aro_labels = parsers.prepare_labels(labels, card)

        le_family = preprocessing.LabelEncoder()
        le_family.fit(amr_family_labels)

        le_aro = preprocessing.LabelEncoder()
        le_aro.fit(aro_labels)

        X_family_train, \
        X_family_test, \
        y_family_train, \
        y_family_test = model_selection.train_test_split(X_family_bitscore_norm,
                                                         amr_family_labels,
                                                         stratify=amr_family_labels,
                                                         test_size=0.33,
                                                         random_state=42)


        X_aro_train, \
        X_aro_test, \
        y_aro_train, \
        y_aro_test = model_selection.train_test_split(X_aro_bitscore_norm,
                                                      aro_labels,
                                                      stratify=aro_labels,
                                                      test_size=0.33,
                                                      random_state=42)

        family_clf = train_family_classifier(X_family_bitscore_norm,
                                             amr_family_labels)

        family_classifiers = train_family_level_classifiers(X_aro_train, y_aro_train, card)


        score(family_clf, family_classifiers, X_family_test, X_aro_test,
              y_family_test, y_aro_test)
       #print(y_family_pred)






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



