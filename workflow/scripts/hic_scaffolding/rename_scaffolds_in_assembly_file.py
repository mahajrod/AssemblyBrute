#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
from RouToolPa.GeneralRoutines import FileRoutines


parser = argparse.ArgumentParser()


parser.add_argument("-a", "--assembly_file", action="store", dest="assembly_file", default=sys.stdin,
                    type=lambda s: FileRoutines.metaopen(s, "r"),
                    help="Input assembly file. Default: stdin.")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    type=lambda s: FileRoutines.metaopen(s, "w"),
                    help="Output file. Default: stdout")
parser.add_argument("--scaffold_syn_file", action="store", dest="scaffold_syn_file", required=True,
                    help="File with scaffold id synonyms")
parser.add_argument("--syn_file_key_column", action="store", dest="syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in synonym file. Default: 0")
parser.add_argument("--syn_file_value_column", action="store", dest="syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in synonym file synonym. Default: 1")

args = parser.parse_args()

syn_df = pd.read_csv(args.scaffold_syn_file, usecols=(args.syn_file_key_column, args.syn_file_value_column),
                     index_col=0, header=None, sep="\t", names=["scaffold_", "syn"])

with FileRoutines.metaopen(args.assembly_file, "r") as in_fd, FileRoutines.metaopen(args.output, "w") as out_fd:
    for line in in_fd:
        if line[0] == ">":
            line_list = line[1:].strip().split()
            line_list[0] = syn_df.loc[line_list[0], "syn"]
            out_fd(">{0}\n".format(" ".join(line_list)))
        else:
            out_fd.write(line)

