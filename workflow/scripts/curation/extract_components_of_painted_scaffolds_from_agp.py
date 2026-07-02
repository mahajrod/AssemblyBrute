#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse

import pandas as pd

from RouToolPa.GeneralRoutines import FileRoutines


parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input painted agp. Default: stdin")

parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix. Required.")

args = parser.parse_args()

agp_df = pd.read_csv(FileRoutines.metaopen(args.input, "r"), sep="\t", header=None,
                     names=["scaffold_id", "start", "end", "part_number", "part_type",
                            "part_id/gap_length", "part_start/gap_type",
                            "part_end/linkage", "orientation/evidence", "comment"],
                     index_col="scaffold_id")

all_contig_series = agp_df[agp_df["part_type"] != "U"]["part_id/gap_length"]
chr_component_series = agp_df[agp_df["comment"] == "Painted"]["part_id/gap_length"]
chr_component_series.to_csv(f"{args.output_prefix}.all_chr.components.ids", sep="\t", header=False, index=False)

for scaffold_id in chr_component_series.index:
    chr_component_series[scaffold_id].to_csv(f"{args.output_prefix}.{scaffold_id}.components.ids", sep="\t", header=False, index=False)
    chr_black_list_series = chr_component_series[~chr_component_series.isin(chr_component_series[scaffold_id])]
    chr_black_list_series.to_csv(f"{args.output_prefix}.{scaffold_id}.pretext.blacklist", sep="\t", header=False, index=False)