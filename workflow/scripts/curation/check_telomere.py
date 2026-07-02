#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
import numpy as np

START_COLUMN_INDEX = 0
END_COLUMN_INDEX = 1
SCORE_COLUMN_INDEX = 2
LENGTH_COLUMN_INDEX = 3
FIVE_PRIME_DIST_COLUMN_INDEX = 4
THREE_PRIME_DIST_COLUMN_INDEX = 5
MIN_END_DISTANCE_COLUMN_INDEX = 6

def get_status(df_row):
    if (df_row.iloc[FIVE_PRIME_DIST_COLUMN_INDEX] >= args.min_distance) and (df_row.iloc[THREE_PRIME_DIST_COLUMN_INDEX] >= args.min_distance):
        return "INTERNAL"
    elif (df_row.iloc[FIVE_PRIME_DIST_COLUMN_INDEX] < args.min_distance) and (df_row.iloc[THREE_PRIME_DIST_COLUMN_INDEX] >= args.min_distance):
        return "FIVE_PRIME"
    elif (df_row.iloc[THREE_PRIME_DIST_COLUMN_INDEX] < args.min_distance) and (df_row.iloc[FIVE_PRIME_DIST_COLUMN_INDEX] >= args.min_distance):
        return "THREE_PRIME"
    else:
        return "AMBIGUOUS"


parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input bedgraph file with telomere fraction. Default: stdin")
parser.add_argument("-f", "--fai_file", action="store", dest="fai_file", required=True,
                    help="Input .fai file. Required.")
parser.add_argument("-d", "--min_distance", action="store", dest="min_distance", default=100000, type=int,
                    help="Minimal distance from end for warning. Default: 100'000")
parser.add_argument("-s", "--score_threshold", action="store", dest="score_threshold", default=0.5, type=float,
                    help="Minimum score (fraction of telomere sequences in window) to report window. Default: '0.5'")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix. Required")

args = parser.parse_args()

fai_df = pd.read_csv(args.fai_file, sep="\t", header=None, names=["scaffold", "length"],
                     usecols=[0, 1], index_col="scaffold").sort_values(by=["length", "scaffold"], ascending=(False, True))

telomere_df = pd.read_csv(args.input, sep="\t", header=None, names=["scaffold", "start", "end", "score"],
                          usecols=[0, 1, 2, 3], index_col="scaffold")

telomere_df["scaffold_length"] = fai_df["length"]
telomere_df["five_prime_dist"] = telomere_df["start"]
telomere_df["three_prime_dist"] = telomere_df["scaffold_length"] - telomere_df["end"]
telomere_df["min_end_distance"] = np.minimum(telomere_df["five_prime_dist"], telomere_df["three_prime_dist"])

telomere_df["status"] = telomere_df.apply(get_status, axis=1)

internal_telomere_warning_file = "{0}.bedgraph".format(args.output_prefix)
telomere_window_status_file = "{0}.status".format(args.output_prefix)
telomere_window_status_filtered_file = "{0}.filtered.status".format(args.output_prefix)
telomere_contig_status_file = "{0}.contig.status".format(args.output_prefix)

telomere_df.to_csv(telomere_window_status_file,
                   sep="\t",
                   index=True,
                   header=False)

telomere_filtered_df = telomere_df[telomere_df["score"] >= args.score_threshold]
#print(telomere_filtered_df)
telomere_filtered_df[telomere_filtered_df["status"] == "INTERNAL"][["start", "end", "score"]].to_csv(internal_telomere_warning_file,
                                                                                                     sep="\t",
                                                                                                     index=True,
                                                                                                     header=False)

#scaffold_status_df = fai_df[["length"]]
#print(telomere_filtered_df[["status", "score"]].groupby(["scaffold", "status"]).count())