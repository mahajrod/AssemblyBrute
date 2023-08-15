#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input bedgraph file with telomere fraction. Default: stdin")
parser.add_argument("-f", "--fai_file", action="store", dest="fai_file", required=True,
                    help="Input .fai file. Required.")
parser.add_argument("-d", "--min_distance", action="store", dest="min_distance", default=300000, type=int,
                    help="Minimal distance from end for warning. Default: 500'000")
parser.add_argument("-s", "--score_threshold", action="store", dest="score_threshold", default=0.5, type=float,
                    help="Minimum score (fraction of telomere sequences in window) to report window. Default: '0.5'")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output. Default: stdout")

args = parser.parse_args()

fai_df = pd.read_csv(args.fai_file, sep="\t", header=None, names=["scaffold", "length"],
                     usecols=[0, 1], index_col="scaffold").sort_values(by=["length", "scaffold"], ascending=(False, True))

telomere_df = pd.read_csv(args.input, sep="\t", header=None, names=["scaffold", "start", "end", "score"],
                          usecols=[0, 1, 2, 3], index_col="scaffold")

telomere_df["length"] = fai_df["length"]

telomere_df["min_end_distance"] = np.minimum(telomere_df["length"] - telomere_df["start"], telomere_df["end"])
telomere_df[(telomere_df["min_end_distance"] >= args.min_distance) & (telomere_df["score"] >= args.score_threshold)][["start", "end", "score"]].to_csv(args.output,
                                                                                                    sep="\t",
                                                                                                    index=True,
                                                                                                    header=False)
