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

    #generate training data
    dataset, labels = model.generate_training_data(card)
    data = model.prepare_data(dataset, labels, card)

    family_classifier = model.GeneFamilyLevelClassifier(data['train']['family']['X'],
                                                              data['train']['family']['y'],
                                                              card)

    family_classifier.train()

    aro_classifiers = model.SubGeneFamilyModel(data['train']['aro']['X'],
                                               data['train']['aro']['y'],
                                               card)

    aro_classifiers.train()
    model.score(family_classifier.clf, aro_classifiers.family_level_classifiers,
            data['test']['family']['X'], data['test']['family']['y'],
            data['test']['aro']['X'], data['test']['aro']['y'])

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
