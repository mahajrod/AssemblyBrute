#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
from RouToolPa.GeneralRoutines import FileRoutines


parser = argparse.ArgumentParser()


parser.add_argument("-c", "--contig_bed", action="store", dest="contig_bed", required=True,
                    help="Bed/bedgraph file with annotations. Required")
parser.add_argument("-t", "--transfer_agp", action="store", dest="transfer_agp", required=True,
                    help="Transfer agp file. Required")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="File to write liftovered file. Default: stdout")

args = parser.parse_args()

transfer_agp_df = pd.read_csv(args.transfer_agp, header=None, sep="\t",
                              names=["assembly", "assembly_start", "assembly_end", "component_number", "component_type",
                                     "contig", "contig_start", "contig_end", "contig_strand"])
transfer_agp_df["assembly_start"] -= 1
transfer_agp_df["contig_start"] -= 1

with FileRoutines.metaopen(args.contig_bed, "r") as in_fd, FileRoutines.metaopen(args.output, "r") as out_fd:
    for line in in_fd:
        if line[0] == "#":
            out_fd.write(line)
        else:
            line_list = line.strip().split("\t")
            shift = transfer_agp_df[transfer_agp_df["contig"] == line_list[0]].loc[0, "assembly_start"]
            strand = transfer_agp_df[transfer_agp_df["contig"] == line_list[0]].loc[0, "assembly_strand"]
            contig_length = transfer_agp_df[transfer_agp_df["contig"] == line_list[0]].loc[0, "contig_end"]
            line_list[0] = "assembly"
            if strand == "+":
                line_list[1] = str(int(line_list[1]) + shift)
                line_list[2] = str(int(line_list[2]) + shift)
            else:
                line_list[1] = str(contig_length - int(line_list[1]) + shift)
                line_list[2] = str(contig_length - int(line_list[2]) + shift)
                line_list[1], line_list[2] = line_list[2], line_list[1]

            out_fd.write("\t".join(line_list) + "\n")
