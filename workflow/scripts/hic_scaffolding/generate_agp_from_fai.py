#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
from RouToolPa.GeneralRoutines import FileRoutines


parser = argparse.ArgumentParser()


parser.add_argument("-f", "--fai_file", action="store", dest="fai_file", default=sys.stdin,
                    type=lambda s: FileRoutines.metaopen(s, "r"),
                    help="Input fai file.")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    type=lambda s: FileRoutines.metaopen(s, "w"),
                    help="Prefix of output files")

args = parser.parse_args()

fai_df = pd.read_csv(args.fai_file, sep="\t", header=None, names=["scaffold", "length"],
                     usecols=[0, 1]).sort_values(by=["length", "scaffold"], ascending=(False, True))
fai_df["scaffold_start"] = 1
fai_df["scaffold_end"] = fai_df["length"]
fai_df["component_number"] = 1
fai_df["component_type"] = "W"
fai_df["strand"] = "+"
fai_df["object"] = fai_df["scaffold"]
fai_df["object_start"] = fai_df["scaffold_start"]
fai_df["object_end"] = fai_df["scaffold_end"]

fai_df[["object", "object_start", "object_end", "component_number", "component_type", "scaffold",
        "scaffold_start", "scaffold_end", "strand"]].to_csv(args.output, sep="\t", header=False, index=False)
