#!/usr/bin/env python
import os

from amrtime import utils
from amrtime import parsers
#from amrtime import database
from amrtime import model

# how do I make params that only evaluate once again?
RANDOM_STATE = 42

def train(args):
    card = parsers.CARD(args.card_fp)

    #generate training data from card.json into labelled fastq
    family_fastq, family_labels, aro_fastq, aro_labels = \
            model.generate_training_data(card, args.redo)

    #encode this prepared data
    #for the family-level model
    family_data = model.prepare_data(family_fastq, family_labels, "family" ,
                                     card, 'training_data/family_data.pkl')

    #and for the aro/subfamily level models
    aro_data = model.prepare_data(subfamily_fastq, subfamily_labels, "aro",
                                 card, 'training_data/subfamily_data.pkl')

    # train the family level classifier
    family_classifier = model.GeneFamilyLevelClassifier(family_data['train']['X'],
                                                        family_data['train']['y'],
                                                              card)
    family_classifier.train()

    # and the aro level classifier
    aro_classifiers = model.SubGeneFamilyModel(aro_data['train']['X'],
                                               aro_data['train']['y'],
                                               card)
    aro_classifiers.train()

    # then evaluate both models
    model.score(family_classifier.clf, aro_classifiers.family_level_classifiers,
            family_data['test']['X'], family_data['test']['y'],
            aro_data['test']['X'], aro_data['test']['y'])

def classify(args):
    card = parsers.CARD(args.card_fp)
    trained_model = parsers.CARD(args.model)




if __name__ == '__main__':

    parser = utils.get_parser()
    args = parser.parse_args()

    # ugh having an issue where the same sequence is appearing multiple times
    # in the card dna sequences thus leading to a mismatch in sequence number
    # when I key on accession during encoding/
    # why are the same sequences in card.dna_sequences multiple times?
    # have I fucked up my code or is that a thing in the CARD database
    if args.mode == 'train':
        train(args)

    elif args.mode == 'classify':

        #if not os.path.exists(arg.model_dir):
        #    print('Must train first')
        pass
        # if read passes DIAMOND filter
        # encode read to kmers/features (maybe do in bulk for efficiency later)
        # run against trained gene family predictor
        # run against gene family model for highest predictor
        # what is max(y_pred) for read
        # evaluate against that gene family level model in the trained model
        #   dictionary
