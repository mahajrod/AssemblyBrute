#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file. Default: stdout")
parser.add_argument("-g", "--convert_to_gbp", action="store_true", dest="convert_to_gbp", default=False,
                    help="Convert number of bases from bp to Gbp. Conflicts with -m/--convert_to_mbp. Default: False.")
parser.add_argument("-m", "--convert_to_mbp", action="store_true", dest="convert_to_mbp", default=False,
                    help="Convert number of bases from bp to Mbp. Conflicts with -g/--convert_to_gbp. Default: False.")
parser.add_argument("-l", "--label_list", action="store", dest="label_list", type=lambda s: s.split(","),
                    help="Comma separated list of the labels corresponding to the input files. Order must be the same."
                         "Default: not set")
parser.add_argument("input", action="store", nargs="+",
                    help="Space separated list of the input files")

args = parser.parse_args()

label_dict = None
if args.label_list is not None:
    if len(args.label_list) != len(args.input):
        raise ValueError("ERROR!!! Number of labels differs from number of files. Something is missing...")
    label_dict = {filename: label for filename, label in zip(args.input, args.label_list)}

df_list = []
for filename in args.input:
    df_list.append(pd.read_csv(filename, sep="\t", skiprows=1, header=None,
                               names=("reads", label_dict[filename] if label_dict is not None else filename),
                               index_col=0).transpose())
    df_list[-1].index.name = "reads"


final_df = pd.concat(df_list)

if args.convert_to_gbp and args.convert_to_mbp:
    raise ValueError("ERROR!!! Both -g/--convert_to_gbp and -m/--convert_to_gbp options are set. "
                     "Choose only one of them...")
elif args.convert_to_gbp:
    final_df["number_of_bases"] = final_df["number_of_bases"].astype(float) / 1000000000
    final_df.rename(columns={"number_of_bases": "number_of_bases,Gbp"}, inplace=True)
elif args.convert_to_mbp:
    final_df["number_of_bases"] = final_df["number_of_bases"].astype(float) / 1000000
    final_df.rename(columns={"number_of_bases": "number_of_bases,Mbp"}, inplace=True)

final_df.to_csv(args.output, header=True, index=True, sep="\t")
