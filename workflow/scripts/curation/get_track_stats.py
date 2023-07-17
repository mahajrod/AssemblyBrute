#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
from RouToolPa.Parsers.Sequence import CollectionSequence
from RouToolPa.GeneralRoutines import FileRoutines

parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input bedgraph file. Default: stdin")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-n", "--normalize", action="store_true", dest="normalize",
                    help="Normalize by length. Default: False")
parser.add_argument("-t", "--threshold_list", action="store", dest="threshold_list", default=[0.05, 0.1, 0.25, 0.5],
                    type=lambda s: list(map(float, s.split(","))),
                    help="Calculate thresholds (fraction of median relative to the median). Default: 0.5,0.1,0.25,0.5")

args = parser.parse_args()

bedgraph_df = pd.read_csv(args.input, sep="\t", header=None, names=["scaffold", "start", "end", "score"], index_col=0)

if args.normalize:
    bedgraph_df["score"] = bedgraph_df["score"].divide(bedgraph_df["end"] - bedgraph_df["start"])

per_scaffold_stats_df = bedgraph_df["score"].groupby("scaffold").agg(["min", "max", "median", "mean"])
per_scaffold_stats_df.to_csv("{0}.per_scaffold.stat".format(args.output_prefix), sep="\t", header=True, index=True)

stats_df = pd.DataFrame()
stats_df["all"] = bedgraph_df["score"].agg(["min", "max", "median", "mean"])
stats_df = stats_df.transpose()
stats_df.index.name = "scaffold"
stats_df.to_csv("{0}.stat".format(args.output_prefix), sep="\t", header=True, index=True)

threshold_list = [0.0]
for threshold in args.threshold_list:
    threshold_list.append(stats_df.loc["all", "median"] * (1 + threshold))
    threshold_list.append(stats_df.loc["all", "median"] * (1 - threshold))

with open("{0}.thresholds".format(args.output_prefix), "w") as out_fd:
    out_fd.write(",".join(list(map(str, sorted(threshold_list)))))
