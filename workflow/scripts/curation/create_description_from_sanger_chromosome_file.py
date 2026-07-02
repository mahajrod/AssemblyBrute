#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input three column file with sanger chromosomes. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file with description column for submission. Default: stdout")

args = parser.parse_args()


chr_df = pd.read_csv(args.input, sep=',', header=None, names=["scaffold_id", "chromosome_id", "is_chromosome"],)

chr_df["location_description"] = chr_df["is_chromosome"].apply(lambda s: "[location=chromosome]" if s == "yes" else "")
chr_df["chromosome_description"] = chr_df["chromosome_id"].apply(lambda s: "[chromosome={0}]".format(s))
chr_df["description"] = chr_df["location_description"] + " " + chr_df["chromosome_description"]

chr_df[["scaffold_id", "description"]].to_csv(args.output, header=False, sep="\t", index=False)
