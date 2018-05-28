#!/usr/bin/env python

import subprocess
import argparse
import os

# how do I make params that only evaluate once again?
RANDOM_STATE = 42

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument()

    #train/test mode
    #card.json used

class GeneFamilyLevelClassifier():
    def __init__(self, card, card):
        self.card = card

    def simulate_reads(self):
        # wait I just called art directly for this when I did it on
        # veles so just get that code when you access it
        # art_illumina

    def rebalancing(self):
        """
        SMOTE with Tomeks link cleaning using imbalanced learning library
        """
        pass

    def build_X_y(self):
        pass

    def train(self):
        #self.X, self.y = self.build_X_y()


class SubGeneFamilyModel(GeneFamilyLevelClassifier):
    """
    Classifier per gene family i.e. attempting to determine a function
    to map from an arbitrary AMR gene family to the specific set of
    AROs within that gene family
    """
    def __init__(self, gene_family, family_to_aro_map, classifier):
        if classifier == 'NaiveBayes':
            from sklearn import naive_bayes
            # only if features are integers i.e. k-mer encoding maybe not
            # some of the other encodings... not that k-mers are independent
            # so violate the assumptions of Naive Bayes
            self.clf = naive_bayes.MultinomialNB()
        else:
            raise NotImplementedError(classifier)




class CARD():
    def __init__(self, card_json_fp):
        with open(card_json_fp) as fh:
            self.card = json.load(card_json_fp)
            self.version = self.card['_version']

            # to avoid having to except them later when parsing the other
            # entries
            del self.card['_version']
            del self.card['_timestamp']
            del self.card['_comment']

            self.aro_to_gene_family = build_aro_to_gene_family()
            self.gene_family_to_aro = build_gene_family_to_aro()

    def build_aro_to_gene_family(self):
        aro_to_gene_family = {}
        for card_item in self.card.values():
            ARO_acc = card_item['ARO_accession']

            # as multiple gene families are possible per ARO
            # although they are relatively rare so maybe I should talk with
            # Andrew about making them unique
            gene_families = []
            for category in card_item['ARO_category'].values():
                if category['category_aro_class_name'] == 'AMR Gene Family':
                    gene_families.append(category['category_aro_name'])
            aro_to_gene_family.update({ARO_acc: gene_family})
        return aro_to_gene_family

    def build_gene_family_to_aro(self):
        gene_family_to_aro = {}
        for key, value in self.aro_to_gene_family.items():
            for gene_family in value:
                if gene_family not in gene_family_to_aro:
                    gene_family_to_aro.update({gene_family: [key]})
                else:
                    gene_family_to_aro[gene_family].append(key)
        return gene_family_to_aro


class Read():
    """
    Class to keep the reads encoding nice and tidy
    """
    def __init__(self, read_id, read_seq):
        self.id = read_id
        self.seq = read_seq

    @param
    def encoded_seq(self):
        self.encoded_seq = encoding(self.seq)

    def encoding(self, sequence):
        """
        k-mer decomposition might be simplest although might have to be
        done at the file level in order to represent all k-mers
        """

if __name__ == '__main__':

    #train mode
        # get gene family to aro map and reverse
        # generate_training for card sequences (can add prev later)
        # parse reads using Read class and labels using split...
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

    #evaluate mode
        # if read passes DIAMOND filter
        # encode read to kmers/features (maybe do in bulk for efficiency later)
        # run against trained gene family predictor
        # run against gene family model for highest predictor
        # what is max(y_pred) for read
        # evaluate against that gene family level model in the trained model
        #   dictionary


