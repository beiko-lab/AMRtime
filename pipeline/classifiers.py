#!/usr/bin/env python


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


