#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
from reprlib import aRepr

import pandas as pd
from RouToolPa.GeneralRoutines import FileRoutines


parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input bed-like file. Default: stdin")
parser.add_argument("-d", "--id_list", action="store", dest="id_list", default=None,
                    help="While with ids to filter. Default: not set, i.e. input will be just forwarded to the output")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="blacklist",
                    help="Mode of the filtration. Allowed: blacklist(default), whitelist")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="File to write output. Default: stdout")

args = parser.parse_args()

if args.id_list is None:
    id_list = None
else:
    try:
        id_list = list(pd.read_csv(args.id_list, sep="\t", header=None).squeeze())
    except pd.errors.EmptyDataError:
        sys.stderr.write("EMPTY ID FILE\n")
        id_list = None
if args.mode == "blacklist":
    comparator_func = lambda s: s not in id_list
elif args.mode == "whitelist":
    comparator_func = lambda s: s in id_list
else:
    raise ValueError(f"ERROR!!! Unknown mode: {args.mode}")

with FileRoutines.metaopen(args.input, "r") as in_fd, FileRoutines.metaopen(args.output, "w") as out_fd:
    if id_list is None:
        for line in in_fd:
            out_fd.write(line)
    else:
        for line in in_fd:
            if comparator_func(line.split("\t")[0]):
                out_fd.write(line)
