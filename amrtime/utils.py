#!/usr/bin/env python
import os, sys
import argparse
import subprocess
import pandas as pd
import shutil
import logging
from pathlib import Path

def check_file(path: str) -> Path:
    """
    Check an input file exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_file():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} can't be read")


def check_dependencies():
    """
    Check all dependencies exist and work
    """
    missing=False
    for program in ["diamond version", "vsearch -v"]:
        try:
            output = subprocess.run(program, shell=True, check=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, encoding='utf-8')
            version = output.stdout
            logging.debug(f"Tool {program.split()[0]} is installed: {version}")
        except:
            logging.error(f"Tool {program.split()[0]} is not installed")
            missing=True
    if missing:
        logging.error("One or more dependencies are missing please install")
        sys.exit(1)
    else:
        logging.info("All dependencies found")


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

