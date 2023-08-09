#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
from RouToolPa.GeneralRoutines import FileRoutines


parser = argparse.ArgumentParser()


parser.add_argument("-a", "--assembly_agp", action="store", dest="assembly_agp", required=True,
                    help="Assembly agp file. Required")
parser.add_argument("-l", "--liftover_agp", action="store", dest="liftover_agp", required=True,
                    help="Liftover agp file. Required")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="File to write transfer agp file. Default: stdout")

args = parser.parse_args()

assembly_agp_df = pd.read_csv(args.assembly_agp, header=None, sep="\t",
                              names=["assembly", "assembly_start", "assembly_end", "component_number", "component_type",
                                     "contig", "contig_start", "contig_end", "contig_strand"])
liftover_agp_df = pd.read_csv(args.liftover_agp, header=None, sep="\t",
                              names=["contig", "contig_start", "contig_end", "component_number", "component_type",
                                     "scaffold", "scaffold_start", "scaffold_end", "strand"])

assembly_agp_df["contig"] = liftover_agp_df["contig"]
assembly_agp_df["contig_start"] = liftover_agp_df["contig_start"]
assembly_agp_df["contig_end"] = liftover_agp_df["contig_end"]

assembly_agp_df.to_csv(args.output, sep="\t", header=False, index=False)
