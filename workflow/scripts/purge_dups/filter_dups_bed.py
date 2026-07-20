#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import sys
import argparse
import pandas as pd

def read_str_series(s):
    if s is None:
        return pd.Series()
    if os.path.exists(s):
        return pd.read_csv(s, header=None, dtype=str).squeeze("columns")
    else:
        return pd.Series(s.split(","))

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Bed file with purge_dups output. Default: stdin")
parser.add_argument("-b", "--type_blacklist", action="store", dest="type_blacklist", default=None,
                    help="Comma-separated list of purge_dups scaffold types to be removed from input dups.bed. "
                         "Default: not set")
parser.add_argument("-w", "--type_whitelist", action="store", dest="type_whitelist", default=None,
                    help="Comma-separated list of purge_dups scaffold types to be retained in the input dups.bed. "
                         "Default: not set")
parser.add_argument("--scaffold_blacklist", action="store", dest="scaffold_blacklist_sr", default=pd.Series(),
                    type=read_str_series,
                    help="Comma-separated list or file with ids of scaffolds to be removed from input dups.bed. "
                         "Scaffold-level filters are applied AFTER type filters. "
                         "Default: not set")
parser.add_argument("--scaffold_whitelist", action="store", dest="scaffold_whitelist_sr", default=pd.Series(),
                    type=read_str_series,
                    help="Comma-separated list or file with ids of scaffolds allowed to be retained in the input dups.bed. "
                         "Scaffold-level filters are applied AFTER type filters. "
                         "Default: not set")

parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Prefix of output files")

args = parser.parse_args()

if args.type_blacklist is not None:
    args.type_blacklist = map(lambda s: s.upper(), args.type_blacklist.split(","))
if args.type_whitelist is not None:
    args.type_whitelist = map(lambda s: s.upper(), args.type_whitelist.split(","))

dups_bed_df = pd.read_csv(args.input, sep="\t", header=None, index_col=0,
                          names=["scaffold", "start", "stop", "type", "duplicate_id"])

print("Filtering purge_dups bed by type...")
if (args.type_blacklist is not None) and (args.type_whitelist is not None):
    print("\tBoth whitelist and blacklist are set for filtering by type...")
    allowed_set = set(args.type_whitelist) - set(args.type_blacklist)
    filtered_df = dups_bed_df[dups_bed_df["type"].isin(allowed_set)]
elif args.type_whitelist is not None:
    print("\tOnly whitelist is set for filtering by type...")
    allowed_set = set(args.type_whitelist)
    filtered_df = dups_bed_df[dups_bed_df["type"].isin(allowed_set)]
elif args.type_blacklist is not None:
    print("\tOnly blacklist is set for filtering by type...")
    blocked_set = set(args.type_blacklist)
    filtered_df = dups_bed_df[~dups_bed_df["type"].isin(blocked_set)]
else:
    print("\tNeither blacklist nor whitelist is set for filtering by type. Skipping...")
    filtered_df = dups_bed_df

print("Filtering purge_dups bed by scaffolds...")
if (not args.scaffold_whitelist_sr.empty()) and (not args.scaffold_blacklist_sr.empty()):
    print("\tBoth whitelist and blacklist are set for filtering by scaffolds...")
    allowed_set = set(args.scaffold_whitelist_sr) - set(args.scaffold_blacklist_sr)
    filtered_df = filtered_df[filtered_df.index.isin(allowed_set)]
elif not args.scaffold_whitelist_sr.empty():
    print("\tOnly whitelist is set for filtering by scaffolds...")
    allowed_set = set(args.scaffold_whitelist_sr)
    filtered_df = filtered_df[filtered_df.index.isin(allowed_set)]
elif not args.scaffold_blacklist_sr.empty():
    print("\tOnly blacklist is set for filtering by scaffolds...")
    blocked_set = set(args.scaffold_blacklist_sr)
    filtered_df = dups_bed_df[~filtered_df.index.isin(blocked_set)]
else:
    print("\tNeither blacklist nor whitelist is set for filtering by scaffolds. Skipping...")

filtered_df.to_csv(args.output,  sep="\t", header=False, index=True)
