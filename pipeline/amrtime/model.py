#!/usr/bin/env python

import os
import subprocess
import numpy as np
from tqdm import tqdm

from imblearn.over_sampling import SMOTE
from sklearn import model_selection
from sklearn import ensemble
from sklearn import preprocessing
from sklearn import dummy
from sklearn import metrics
from sklearn.externals import joblib
import pickle

from amrtime import encoding
from amrtime import parsers
from amrtime import utils

class GeneFamilyLevelClassifier():
    def __init__(self, X_family_train, y_family_train, card, name=None):
        self.card = card
        self.X = X_family_train
        self.y = y_family_train
        if name is None:
            self.name = "model"
        else:
            self.name = name

    def train(self):
        print("Training Family Classifier")
        print("Rebalancing Data")
        if os.path.exists('models/family.pkl'):
            self.clf = joblib.load('models/family.pkl')
        else:
            X_resampled, y_resampled = SMOTE(kind='borderline1').fit_sample(self.X,
                                                                            self.y)
            print("Training Model")
            clf = ensemble.RandomForestClassifier()
            clf.fit(X_resampled, y_resampled)
            joblib.dump(clf, 'models/family.pkl')
            self.clf = clf

    def test(self):
        print("Evaluating


class SubGeneFamilyModel(GeneFamilyLevelClassifier):
    """
    Classifier per gene family i.e. attempting to determine a function
    to map from an arbitrary AMR gene family to the specific set of
    AROs within that gene family
    """
    def __init__(self, X_aro_train, y_aro_train, card):
        self.X = X_aro_train
        self.y = y_aro_train
        self.card = card

    def train(self):
        """
        Train the family level classifiers
        """
        print("Training ARO classifiers for each family")
        family_level_classifiers = {}
        for family in tqdm(self.card.gene_family_to_aro.keys()):

            family_name = family.replace(' ', '_').replace('/', '_')
            # get all the aros relevant to the family
            family_aros = self.card.gene_family_to_aro[family]

            # filter input to just the columns containing similarity to aros
            # within the family
            X_train = self.X[family_aros]

            # get the indices where the label is one of the AROs belonging to
            # the current family
            label_indices = [ix for ix, x in enumerate(self.y) if x in family_aros]

            y_train = np.array(self.y)[label_indices]

            # grab only the reads where the label index is an ARO belonging to the
            # current family being trained
            X_train = X_train.iloc[label_indices]

            if os.path.exists('models/{}.pkl'.format(family_name)):
                family_clf = joblib.load('models/{}.pkl'.format(family_name))
                family_level_classifiers.update({family: [family_clf, family_aros]})
                continue

            # i.e. if family only has a single member
            if len(family_aros) == 1:
                family_clf = dummy.DummyClassifier(strategy='constant', constant=family_aros[0])
                family_clf.fit(X_train, y_train)

                joblib.dump(family_clf, 'models/{}.pkl'.format(family_name))
                family_level_classifiers.update({family: [family_clf, family_aros]})
            else:
                # rebalance using SMOTE
                X_resampled, y_resampled = SMOTE(kind='borderline1').fit_sample(X_train, y_train)

                family_clf = ensemble.RandomForestClassifier()

                family_clf.fit(X_resampled, y_resampled)

                joblib.dump(family_clf, 'models/{}.pkl'.format(family_name))
                family_level_classifiers.update({family: [family_clf, family_aros]})

        self.family_level_classifiers = family_level_classifiers


def generate_training_data(card):
    """
    Simulate training data from the card.json data
    """

    if not os.path.exists('training_data'):
        os.mkdir('training_data')

    if not os.path.exists('training_data/card_proteins.faa'):
        card.write_proteins('training_data/card_proteins.faa')

    if not os.path.exists('training_data/card_nucleotides.fna'):
        card.write_nucleoties('training_data/card_nucleotides.fna')

    # generate training data
    if not os.path.exists('training_data/metagenome.fq'):
        subprocess.check_call('art_illumina -q -na -ef -sam -ss MSv3 -i training_data/card_nucleotides.fna -f 10 -l 250 -rs 42 -o training_data/metagenome', shell=True)
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


def split_data(X_aro, aro_labels, aro_encoder,
               X_family, amr_family_labels, family_encoder):
    """
    Stratified split on AROs
    """

    indices = np.arange(0, len(aro_labels))

    X_aro_train, \
    X_aro_test, \
    y_aro_train, \
    y_aro_test, \
    indices_train, \
    indices_test = model_selection.train_test_split(X_aro,
                                                    aro_labels,
                                                    indices,
                                                    stratify=aro_labels,
                                                    test_size=0.20,
                                                    random_state=42)

    X_family_train = X_family.iloc[indices_train]
    X_family_test = X_family.iloc[indices_test]
    y_family_train = [amr_family_labels[x] for x in indices_train]
    y_family_test = [amr_family_labels[x] for x in indices_test]


    data = {'train': {'family': {'X': X_family_train,
                                 'y': y_family_train},
                      'aro': {'X': X_aro_train,
                              'y': y_aro_train}},
            'test': {'family': {'X': X_family_test,
                                 'y': y_family_test},
                      'aro': {'X': X_aro_test,
                              'y': y_aro_test}},
            'encoders' : {'family': family_encoder,
                           'aro': aro_encoder}}
    return data



def score(family_clf, family_classifiers, X_family, y_family, X_aro, y_aro):

        # predict and calculate performance for X
        print('Predicting Families')
        y_family_pred = family_clf.predict(X_family)

        utils.classification_report_csv(metrics.classification_report(y_family,
                                                                y_family_pred),
                                  "family_test_report.tsv")

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

        utils.classification_report_csv(metrics.classification_report(np.array(y_aro)[all_indices],
                                                                y_aro_pred),
                                  "aro_test_report.tsv")

def prepare_data(dataset, labels, card):
    # run diamond filter to get homology encoding
    homology_encoding = encoding.Homology(dataset,
            'training_data/card_proteins.faa', 'DIAMOND')
    X_family, X_aro = homology_encoding.encode(card, 'bitscore', norm=True)

    #X_family_bitscore_norm = pd.read_pickle('family_bitscore.pkl')
    #X_aro_bitscore_norm = pd.read_pickle('aro_bitscore.pkl')

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

    data = split_data(X_aro, aro_labels, le_aro, X_family, amr_family_labels,
                     le_family)

    joblib.dump(data, 'training_data/data.pkl')

    return data


