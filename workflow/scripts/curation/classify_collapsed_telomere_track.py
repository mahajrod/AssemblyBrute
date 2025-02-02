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
MEDIAN_COLUMN_INDEX = 4
MEANCOLUMN_INDEX = 5
STDEV_COLUMN_INDEX = 6
MODE_COLUMN_INDEX = 7
ABSMIN_COLUMN_INDEX = 8
FIVE_PRIME_DIST_COLUMN_INDEX = 9
THREE_PRIME_DIST_COLUMN_INDEX = 10

def get_status(df_row):
    if (df_row.iloc[FIVE_PRIME_DIST_COLUMN_INDEX] >= args.max_distance) and (df_row.iloc[THREE_PRIME_DIST_COLUMN_INDEX] >= args.max_distance):
        return "INTERNAL"
    elif (df_row.iloc[FIVE_PRIME_DIST_COLUMN_INDEX] < args.max_distance) and (df_row.iloc[THREE_PRIME_DIST_COLUMN_INDEX] >= args.max_distance):
        return "FIVE_PRIME"
    elif (df_row.iloc[THREE_PRIME_DIST_COLUMN_INDEX] < args.max_distance) and (df_row.iloc[FIVE_PRIME_DIST_COLUMN_INDEX] >= args.max_distance):
        return "THREE_PRIME"
    else:
        return "AMBIGUOUS"


parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input bedgraph file with telomere fraction. Default: stdin")
parser.add_argument("-f", "--fai_file", action="store", dest="fai_file", required=True,
                    help="Input .fai file. Required.")
parser.add_argument("-d", "--max_distance", action="store", dest="max_distance", default=10000, type=int,
                    help="Minimal distance from end for warning. Default: 10'000")
parser.add_argument("-s", "--score_threshold", action="store", dest="score_threshold", default=0.6, type=float,
                    help="Minimum score of candidate telomeric region Default: '0.6'")
parser.add_argument("-t", "--score_type", action="store", dest="score_type", default="median",
                    help="Type of the score to use for filtering. Allowed: 'mean', 'median'(default), 'mode', 'absmean'")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix. Required")

args = parser.parse_args()

fai_df = pd.read_csv(args.fai_file, sep="\t", header=None, names=["scaffold", "length"],
                     usecols=[0, 1], index_col="scaffold").sort_values(by=["length", "scaffold"], ascending=(False, True))

telomere_df = pd.read_csv(args.input, sep="\t", header=0,
                          usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8], index_col="#scaffold")
#telomere_df.index.name = "scaffold"

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
                   header=True)

telomere_filtered_df = telomere_df[telomere_df[args.score_type] >= args.score_threshold]
#scaffold_status_df = fai_df[["length"]]
print(telomere_filtered_df[["status", args.score_type]].groupby(["#scaffold", "status"]).count())
print(telomere_filtered_df[["status", args.score_type]])
#print(telomere_filtered_df)
#telomere_filtered_df[telomere_filtered_df["status"] == "INTERNAL"][["start", "end", "score"]].to_csv(internal_telomere_warning_file,
#                                                                                                     sep="\t",
#                                                                                                     index=True,
#                                                                                                     header=False)

