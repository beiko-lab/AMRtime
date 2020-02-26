#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.1.0"

import time
import logging



from amrtime import utils
from amrtime import parsers
#from amrtime import database
from amrtime import model

# how do I make params that only evaluate once again?
RANDOM_STATE = 42


def check_dependencies():
    """
    Check all dependencies exist and work
    """

    missing=False
    for program in ["diamond version", "vsearch -v"]:
        try:
            output = subprocess.run(program, shell=True, check=True,
                                    stdout=subprocess.PIPE, encoding='utf-8')
            version = output.stdout.split('\n')[0]
            logging.debug(f"Tool {program.split()[0]} is installed {version}")
        except:
            logging.error(f"Tool {program.split()[0]} is not installed")
            missing=True
    if missing:
        logging.error("One or more dependencies are missing please install")
        sys.exit(1)
    else:
        logging.debug("All dependencies found")


def run(args):
    """
    Main runner function for AMRtime

    Parameters:
        args: arguments parsed from argparse
    """
    if args.run-mode == 'train':
        run_name = 'Train'
    elif args.run-mode == 'predict':
        run_name = 'Predict'

    if args.verbose:
        logging.basicConfig(format="%(levelname)s:%(message)s",
                            level=logging.DEBUG,
                            handlers=[logging.FileHandler(f"{run_name}.log"),
                                logging.StreamHandler()])

    else:
        logging.basicConfig(format="%(levelname)s:%(message)s",
                            level=logging.INFO,
                            handlers=[logging.FileHandler(f"{run_name}.log"),
                                logging.StreamHandler()])

    #logging.info(f"Started AMRtime
    check_dependencies()


def train(args):
    """
    Train the amrtime models from a supplied CARD released using the
    card.json file
    """

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
