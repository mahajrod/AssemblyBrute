#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Bed file with purge_dups output. Default: stdin")
parser.add_argument("-b", "--blacklist", action="store", dest="blacklist", default=None,
                    help="Comma-separated list of purge_dups scaffold types to be "
                         "removed from input dups.bed. Default: not set")
parser.add_argument("-w", "--whitelist", action="store", dest="whitelist", default=None,
                    help="Comma-separated list of purge_dups scaffold types to be "
                         "retained in the input dups.bed. Default: not set")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Prefix of output files")

args = parser.parse_args()

if args.blacklist is not None:
    args.blacklist = map(lambda s: s.upper(), args.blacklist.split(","))
if args.whitelist is not None:
    args.whitelist = map(lambda s: s.upper(), args.whitelist.split(","))

dups_bed_df = pd.read_csv(args.input, sep="\t", header=None, index_col=0,
                          names=["scaffold", "start", "stop", "type", "duplicate_id"])
if (args.blacklist is not None) and (args.whitelist is not None):
    allowed_set = set(args.whitelist) - set(args.blacklist)
    filtered_df = dups_bed_df[dups_bed_df["type"].isin(allowed_set)]
elif args.whitelist is not None:
    allowed_set = set(args.whitelist)
    filtered_df = dups_bed_df[dups_bed_df["type"].isin(allowed_set)]
elif args.blacklist is not None:
    blocked_set = set(args.blacklist)
    filtered_df = dups_bed_df[~dups_bed_df["type"].isin(blocked_set)]
else:
    dups_bed_df.to_csv(args.output,  sep="\t", header=False, index=True)
    exit(0)
filtered_df.to_csv(args.output,  sep="\t", header=False, index=True)
