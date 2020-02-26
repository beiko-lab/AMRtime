#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse

from amrtime import amrtime
from amrtime import utils

if __name__ == '__main__':

    parent_parser = argparse.ArgumentParser(description="AMRtime metagenomic AMR detector",
            prog="AMRtime")
    parent_parser.add_argument('-v', '--version', action='version',
                        version=f"%(prog)s {amrtime.__version__}",
                        help='Display version and exit')
    parent_parser.add_argument('--verbose', action='store_true', default=False,
                        help="Run with verbose output")
    parent_parser.add_argument('-j', '--num_threads', default=1, type=int,
                        help="Number of processor threads to use")


    subparsers = parent_parser.add_subparsers(help='AMRtime run-mode',
                                              dest="run-mode")
    subparsers.required = True

    parser_train = subparsers.add_parser('train', help='Train AMRtime models',
                                        parents=[parent_parser],
                                        add_help=False)
    parser_train.add_argument("-c", "--card", required=True,
                        type=lambda x: utils.is_valid_file(parser, x),
                        help="Path to CARD json file")
    parser_train.add_argument("-o", "--output",
                              help="Model generation output folder")

    #parser_train.add_argument("--redo", action="store_true",
    #                    help="force overwrite and regeneration of data")

    parser_predict = subparsers.add_parser('predict',
                                    help='Predict AMR genes from a metagenome',
                                    parents=[parent_parser],
                                    add_help=False)

    parser_predict.add_argument("-m", "--model_dir", required=True,
                                help="Folder containing trained model")
    parser_predict.add_argument("-i", "--input", required=True,
                                type=lambda x: utils.is_valid_file(parser, x),
                                nargs="+",
                                help="Input metagenome in fastq format")
    parser_predict.add_argument("-o", '--output',
                                help="Prediction output folder")

    args = parent_parser.parse_args()

    #amrtime.run(args)
