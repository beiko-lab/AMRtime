#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.1.0"

import time
import logging
import sys
import subprocess
import os

from amrtime import utils
from amrtime import parsers
#from amrtime import database
from amrtime import model

# how do I make params that only evaluate once again?
RANDOM_STATE = 42


def run(args):
    """
    Main runner function for AMRtime
    Sets up logging checks dependencies then calls the appropriate
    training or prediction function

    Parameters:
        args: arguments parsed from argparse
    """
    # generating a run name for logging purposes and output purposes
    if args.output_folder:
        run_name = args.mode + " " + os.path.basename(args.output_folder)
    else:
        run_name = args.mode + " " + str(int(time.time()))
        args.output_folder = run_name.replace(' ', '_')


    # check output folder specified
    if os.path.exists(args.output_folder) and not args.force:
        print(f"Output folder: '{args.output_folder}' already exists "
               "please remove or use --force to overwrite before running",
               file=sys.stderr)
        sys.exit(1)
    elif os.path.exists(args.output_folder) and args.force:
        utils.clean_and_remake_folder(args.output_folder)
    else:
        os.mkdir(args.output_folder)


    if args.verbose:
        logging.basicConfig(format="%(levelname)s:%(message)s",
                            level=logging.DEBUG,
                            handlers=[logging.FileHandler(f"{args.output_folder}/{run_name}.log"),
                                logging.StreamHandler()])

    else:
        logging.basicConfig(format="%(levelname)s:%(message)s",
                            level=logging.INFO,
                            handlers=[logging.FileHandler(f"{args.output_folder}/{run_name}.log"),
                                logging.StreamHandler()])

    logging.info(f"Started AMRtime '{run_name}' outputting to {args.output_folder}")
    utils.check_dependencies()

    if args.mode == 'train':
        train(args)
    elif args.mode =='predict':
        predict(args)


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


def predict(args):
    pass

    # check model_dir for trained models and card metadata
    #card = parsers.CARD(args.card_fp)
    #trained_model = parsers.CARD(args.model)

    # if read passes DIAMOND filter
    # encode read to kmers/features (maybe do in bulk for efficiency later)
    # run against trained gene family predictor
    # run against gene family model for highest predictor
    # what is max(y_pred) for read
    # evaluate against that gene family level model in the trained model
    #   dictionary

    # ugh having an issue where the same sequence is appearing multiple times
    # in the card dna sequences thus leading to a mismatch in sequence number
    # when I key on accession during encoding/
    # why are the same sequences in card.dna_sequences multiple times?
    # have I fucked up my code or is that a thing in the CARD database
        #if not os.path.exists(arg.model_dir):
        #    print('Must train first')

