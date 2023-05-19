#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse

from amrtime import amrtime
from amrtime import utils

def main():
    """
    Console script for amrtime
    """
    parent_parser = argparse.ArgumentParser(description="AMRtime metagenomic AMR detector",
                                            prog="amrtime")

    parent_parser.add_argument('-v', '--version', action='version',
                        version=f"%(prog)s {amrtime.__version__}",
                        help='Display version and exit')
    parent_parser.add_argument('--verbose', action='store_true', default=False,
                        help="Run with verbose output")
    parent_parser.add_argument('-j', '--num_threads', default=1, type=int,
                        help="Number of processor threads to use")
    parent_parser.add_argument('--force', default=False, action='store_true',
                        help="Force overwrite of existing output folder")

    subparsers = parent_parser.add_subparsers(help='amrtime run-mode',
                                              dest="mode",
                                              required=True)

    parser_train = subparsers.add_parser('train',
            help='Train amrtime models from a CARD release')
    parser_train.add_argument("-o", "--output_folder",
                              help="Folder to output trained model")
    parser_train.add_argument("-c", "--card", required=True,
                        type=utils.check_file,
                        help="Path to CARD.json file")

    parser_predict = subparsers.add_parser('predict',
                                    help='Predict AMR genes in a metagenome')
    parser_predict.add_argument("-o", '--output_folder',
                                help="Prediction output folder")
    parser_predict.add_argument("-m", "--model_dir", required=True,
                                help="Folder containing trained model")
    parser_predict.add_argument("-i", "--input", required=True,
                                type=utils.check_file,
                                nargs="+",
                                help="Input metagenome in fastq format")

    args = parent_parser.parse_args()

    amrtime.run(args)

if __name__ == '__main__':
    sys.exit(main())
