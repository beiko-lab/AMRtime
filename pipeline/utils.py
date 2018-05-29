#!/usr/bin/env python

import os
import argparse

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file {} does not exist".format(arg))
    else:
        return arg

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=['test', 'train'],
                         help="train model on new card version or evaluate "
                              " trained models on new data")
    parser.add_argument("card_fp", type=lambda x: is_valid_file(parser, x),
                         help="Path to CARD json file")

    return parser

