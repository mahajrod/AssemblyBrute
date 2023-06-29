#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
from RouToolPa.Parsers.Sequence import CollectionSequence
from RouToolPa.GeneralRoutines import FileRoutines

parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input bedgraph file with telomere fraction. Default: stdin")
parser.add_argument("-f", "--fai_file", action="store", dest="fai_file", required=True,
                    help="Input .fai file. Required.")
parser.add_argument("-d", "--min_distance", action="store", dest="min_distance", default=500000, type=int,
                    help="Minimal distance from end for warning. Default: 500'000")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output. Default: stdout")

args = parser.parse_args()

fai_df = pd.read_csv(args.fai_file, sep="\t", header=None, names=["scaffold", "length"],
                     usecols=[0, 1], index_col="scaffold").sort_values(by=["length", "scaffold"], ascending=(False, True))

telomere_df = pd.read_csv(args.input, sep="\t", header=None, names=["scaffold", "start", "end", "score"],
                          usecols=[0, 1, 2, 3], index_col="scaffold")
telomere_df["length"] = fai_df["length"]

print(telomere_df)




