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
parser.add_argument("-g", "--transfer_agp",  action="store", dest="transfer_agp", default=None,
                    help="File to write transfer agp. Default: not set")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="File to write output chain file. Default: stdout")

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

assembly_length = assembly_agp_df["assembly_end"].iloc[-1]
assembly_agp_df["assembly_length"] = assembly_length
if args.transfer_agp:
    assembly_agp_df.to_csv(args.transfer_agp, sep="\t", header=False, index=False)

assembly_agp_df["score"] = assembly_agp_df["contig_end"]
assembly_agp_df["chain"] = "chain"
assembly_agp_df["assembly_strand"] = "+"
chain_agp_df = assembly_agp_df[["chain", "score",
                                "assembly", "assembly_length", "assembly_strand", "assembly_start", "assembly_end",
                                "contig", "contig_end", "contig_strand", "contig_start", "contig_end",
                                "component_number"]]
chain_agp_df["assembly_start"] = chain_agp_df["assembly_start"] - 1  # chain file is zero-based in python notation, agp file - is one-based
chain_agp_df["contig_start"] = chain_agp_df["contig_start"] - 1

with FileRoutines.metaopen(args.output, "w") as out_fd:
    for row_tuple in chain_agp_df.itertuples(index=False):
        out_fd.write(" ".join(map(str, row_tuple)) + "\n{0}\n\n".format(row_tuple[8]))


