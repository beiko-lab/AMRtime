#!/usr/bin/env python
import os
import argparse
import subprocess
import pandas as pd
import shutil

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file {} does not exist".format(arg))
    else:
        return arg

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--redo", action="store_true",
                        help="force overwrite and regeneration of data")

    parser.add_argument("mode", choices=['test', 'train'],
                        help="train model on new card version or evaluate "
                              " trained models on new data")

    parser.add_argument("card_fp", type=lambda x: is_valid_file(parser, x),
                        help="Path to CARD json file")
    return parser

def classification_report_csv(report, fp):
    report_data = []
    lines = report.split('\n')
    for line in lines[2:-3]:
        row = {}
        row_data = line.split('      ')
        row['class'] = row_data[-5]
        row['precision'] = float(row_data[-4])
        row['recall'] = float(row_data[-3])
        row['f1_score'] = float(row_data[-2])
        row['support'] = float(row_data[-1])
        report_data.append(row)
    dataframe = pd.DataFrame.from_dict(report_data)
    #print(dataframe)
    dataframe.to_csv(fp, index = False)

def count_reads_per_file(fp):
    line_count = subprocess.check_output(f"wc -l {fp}", shell=True)
    line_count = line_count.decode().split()[0]
    return int(line_count) / 4

def clean_and_remake_folder(fp):
    if os.path.exists(fp):
        shutil.rmtree(fp)
    os.mkdir(fp)

