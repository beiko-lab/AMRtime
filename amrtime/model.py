#!/usr/bin/env python

import os
import shlex
import shutil
import glob
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
    """
    Class containing the classifier for reads into respective AMR families

    - train
    - test
    """
    def __init__(self, card, name=None):
        self.card = card
        if name is None:
            self.name = "model"
        else:
            self.name = name

    def train(self, encoded_reads, family_labels):
        print("Training Family Classifier")
        print("Rebalancing Data")
        if os.path.exists(f'{self.name}/family.pkl'):
            self.clf = joblib.load(f'{self.name}/family.pkl')
        else:
            X_resampled, y_resampled = SMOTE(kind='borderline1').fit_sample(self.X,
                                                                            self.y)
            print("Training Model")
            clf = ensemble.RandomForestClassifier()
            clf.fit(X_resampled, y_resampled)
            joblib.dump(clf, f'{name}/family.pkl')
            self.clf = clf

    def test(self):
        print("Evaluating")



class SubGeneFamilyModel(GeneFamilyLevelClassifier):
    """
    Classifier per gene family i.e. attempting to determine a function
    to map from an arbitrary AMR gene family to the specific set of
    AROs within that gene family
    """
    def __init__(self, card, name=None):
        #self.X = X_aro_train
        #self.y = y_aro_train
        self.card = card
        if name is None:
            self.name = "model"
        else:
            self.name = name


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

            if os.path.exists('{}/{}.pkl'.format(self.name, family_name)):
                family_clf = joblib.load('{}/{}.pkl'.format(self.name, family_name))
                family_level_classifiers.update({family: [family_clf, family_aros]})

            # i.e. if family only has a single member
            if len(family_aros) == 1:
                family_clf = dummy.DummyClassifier(strategy='constant', constant=family_aros[0])
                family_clf.fit(X_train, y_train)

                joblib.dump(family_clf, '{}/{}.pkl'.format(self.name, family_name))
                family_level_classifiers.update({family: [family_clf, family_aros]})
            else:
                # rebalance using SMOTE
                X_resampled, y_resampled = SMOTE(kind='borderline1').fit_sample(X_train, y_train)

                family_clf = ensemble.RandomForestClassifier()

                family_clf.fit(X_resampled, y_resampled)

                joblib.dump(family_clf, f'{self.name}/{family_name}.pkl')
                family_level_classifiers.update({family: [family_clf, family_aros]})

        self.family_level_classifiers = family_level_classifiers


def generate_training_data(card, redo=False):
    """
    Simulate training data from the card.json data

    - sample reads from each family at equal coverage (5x?)
    - read files and determine the number of reads in each family fastq
    - generate additional reads proportional to the missing number
    - resulting in each numbers of reads for each family
    """

    # check if training data has already been generated, if so return the
    # relevant paths unless the redo flag has been used
    if os.path.exists('training_data/family_metagenome.fq') and \
            os.path.exists('training_data/family_labels.tsv') and not redo:
        print("Training data already generated in training_data folder")
        print("Use --redo flag to rebuild")
        return "training_data/metagenome.fq" "training_data/metagenome_labels.tsv"

    if not os.path.exists('training_data'):
        os.mkdir('training_data')


    training_fasta_folder = "training_data/family_seqs"
    if not os.path.exists(training_fasta_folder):
        os.mkdir(training_fasta_folder)

    amr_family_fasta_locs = card.write_nucleotide_families(training_fasta_folder)

    # generate fastq files with even coverage over each family for training
    # subgene level models
    subfamily_fastq_folder = training_fasta_folder + "/subfamily_fq"
    utils.clean_and_remake_folder(subfamily_fastq_folder)

    print("Generating balanced training data for subfamily classifiers")
    subfamily_fastq_locs = {}
    for family_name, family_fp in amr_family_fasta_locs.items():

        fq_fp = card.convert_amr_family_to_filename(family_name,
                          folder=subfamily_fastq_folder)

        with open(os.devnull, 'w') as null_fh:
            subprocess.check_call(f"art_illumina -q -na -ss MSv3 -i {shlex.quote(family_fp)} \
                -f 10 -l 250 -rs 42 -o {shlex.quote(fq_fp)}", shell=True,
                stderr=subprocess.STDOUT, stdout=null_fh)
        subfamily_fastq_locs.update({family_name: fq_fp + ".fq"})


    # cluster the seqs in each family at 99% identity
    print("Clustering sequences within families for family level classifier")
    family_clustered_folder = training_fasta_folder + "/clustered"
    utils.clean_and_remake_folder(family_clustered_folder)
    clustered_fasta_locs = {}
    for family_name, family_fp in amr_family_fasta_locs.items():
        clustered_fp = card.convert_amr_family_to_filename(family_name,
                          folder=family_clustered_folder)

        with open(os.devnull, 'w') as null_fh:
            subprocess.check_call(f"cd-hit -i {shlex.quote(family_fp)} -c 0.99 \
                                    -o {shlex.quote(clustered_fp)}",
                                  shell=True,
                                  stderr=subprocess.STDOUT, stdout=null_fh)
        clustered_fasta_locs.update({family_name: clustered_fp})

        # add suffixes to the clustered sequence filenames so there isn't
        # name collisions
        with open(clustered_fp) as fh:
            seqs = fh.readlines()
        with open(clustered_fp, 'w') as fh:
            for line in seqs:
                if line.startswith(">"):
                    line = line.strip() + "_clust" + '\n'
                    fh.write(line)
                else:
                    fh.write(line)

    # for each of these generate a set of reads at even coverage
    print("Generating training data for family level classifier")
    family_fastq_folder = family_clustered_folder + "/fq"
    utils.clean_and_remake_folder(family_fastq_folder)
    family_read_counts = {}
    family_fastq_locs = {}
    for family_name, clustered_fp in clustered_fasta_locs.items():
        clustered_fq_fp = card.convert_amr_family_to_filename(family_name,
                                            folder=family_fastq_folder)
        with open(os.devnull, 'w') as null_fh:
            subprocess.check_call(f"art_illumina -q -na -ss MSv3 -i {shlex.quote(clustered_fp)} \
                    -f 10 -l 250 -rs 42 -o {shlex.quote(clustered_fq_fp)}", shell=True,
                    stderr=subprocess.STDOUT, stdout=null_fh)

        family_fastq_locs.update({family_name: clustered_fq_fp + ".fq"})
        family_read_counts.update({family_name: utils.count_reads_per_file(\
                    shlex.quote(clustered_fq_fp + ".fq"))})

    # now we find the largest read count
    biggest_count = max([count for count in family_read_counts.values()])

    # for each family generate a new file with the number of reads needed
    # this is currently being done for the unclustered family seqs
    print("Balancing training data for family level classifier")
    family_fastq_topup_locs = {}

    for family_name, family_fp in amr_family_fasta_locs.items():
        top_up_fq_fp = card.convert_amr_family_to_filename(family_name, folder=family_fastq_folder) + "_topup"

        diff_count = biggest_count - family_read_counts[family_name]

        # difference needed per read to get that count as -c generates
        # reads per sequence in fasta
        family_size = len(card.gene_family_to_aro[family_name])
        diff_reads = round(diff_count / family_size)

        if diff_reads > 0:
            with open(os.devnull, 'w') as null_fh:
                try:
                    subprocess.check_call(f'art_illumina -q -na -ss MSv3 -i {shlex.quote(family_fp)} \
                    -c {diff_reads} -l 250 -rs 42 -o {shlex.quote(top_up_fq_fp)}', shell=True,
                        stderr=subprocess.STDOUT, stdout=null_fh)
                # some sequences e.g. dfr have representatives <250 bp long
                # this apparently isn't a problem with fold coverage for unclear
                # reasons buet causes an exception for option of reads per sequence
                # work-around retry with 225 length reads.
                except:
                    try:
                        subprocess.check_call(f'art_illumina -q -na -ss MSv3 -i {shlex.quote(family_fp)} \
                            -c {diff_reads} -l 225 -rs 42 -o {shlex.quote(top_up_fq_fp)}', shell=True,
                            stderr=subprocess.STDOUT, stdout=null_fh)
                    except:
                        subprocess.check_call(f'art_illumina -q -na -ss MSv3 -i {shlex.quote(family_fp)} \
                            -c {diff_reads} -l 190 -rs 42 -o {shlex.quote(top_up_fq_fp)}', shell=True,
                            stderr=subprocess.STDOUT, stdout=null_fh)


            family_fastq_topup_locs.update({family_name: top_up_fq_fp + ".fq"})

    print("Collating training data for family level classifier")
    comb_family_fastq_folder = training_fasta_folder + "/family_fq"
    utils.clean_and_remake_folder(comb_family_fastq_folder)

    family_combined_fastq_locs = {}
    for family_name, family_fastq_fp in family_fastq_locs.items():
        # biggest family shouldn't have a top up file as its already at the
        # maximum size
        if family_name not in family_fastq_topup_locs:
            if family_read_counts[family_name] == biggest_count:
                top_up_fastq_fp = ""
            else:
                print("Bug with read counts for ", family_name)
        else:
            top_up_fastq_fp = family_fastq_topup_locs[family_name]

        combined_fp = card.convert_amr_family_to_filename(family_name,
                                                          folder=comb_family_fastq_folder)

        if len(top_up_fastq_fp) > 0:
            subprocess.check_call(f"cat {shlex.quote(family_fastq_fp)} {shlex.quote(top_up_fastq_fp)} > {shlex.quote(combined_fp)}",
                                  shell=True)
        else:
            subprocess.check_call(f"cat {shlex.quote(family_fastq_fp)} > {shlex.quote(combined_fp)}",
                                  shell=True)

        family_combined_fastq_locs.update({family_name: combined_fp})

    ## combine all family fq into one big file
    paths = " ".join([shlex.quote(fp) for fp in family_combined_fastq_locs.values()])
    family_fastq = "training_data/family_metagenome.fq"
    subprocess.check_call(f"cat {paths} > {shlex.quote(family_fastq)}",
                          shell=True)

    ## build labels
    print("Building Family Labels")
    family_labels = "training_data/family_labels.tsv"
    subprocess.check_call(f"grep '^@gb' {shlex.quote(family_fastq)} > training_data/read_names",
                         shell=True)

    labels_fh = open(family_labels, 'w')
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

    print("Building Subfamily Labels")
    #todo

    return family_fastq, family_labels


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
                                        'training_data/card_proteins.faa',
                                        'DIAMOND')

    X_family, X_aro = homology_encoding.encode(card, 'bitscore', norm=True)

    amr_family_labels, aro_labels = parsers.prepare_labels(labels, card)

    le_family = preprocessing.LabelEncoder()
    le_family.fit(amr_family_labels)

    le_aro = preprocessing.LabelEncoder()
    le_aro.fit(aro_labels)

    data = split_data(X_aro, aro_labels, le_aro, X_family, amr_family_labels,
                     le_family)

    joblib.dump(data, 'training_data/data.pkl')

    return data


