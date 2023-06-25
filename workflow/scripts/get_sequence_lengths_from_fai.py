#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input .fai ")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file - default: stdout")
args = parser.parse_args()

len_df = pd.read_csv(args.input, sep="\t", header=None, names=["scaffold", "length"],
                     index_col="scaffold", usecols=[0, 1]).sort_values(by=["length", "scaffold"], ascending=(False, True))

len_df.to_csv(args.output, sep="\t", header=False, index=True)
