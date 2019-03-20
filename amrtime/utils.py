#!/usr/bin/env python
import os
import argparse
import pandas as pd

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

