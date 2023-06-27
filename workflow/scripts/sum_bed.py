#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
from RouToolPa.GeneralRoutines import FileRoutines


parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input bed file. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output bed file. Default: stdout")
parser.add_argument("-c", "--column", action="store", dest="column", default=3, type=int,
                    help="Column to sum for same regions, 0-based. Default: 3")

args = parser.parse_args()

with FileRoutines.metaopen(args.input, "r") as in_fd, FileRoutines.metaopen(args.output, "w") as out_fd:
    prev_line_list = in_fd.readline().strip().split("\t")
    prev_line_list = prev_line_list[0:3] + [float(prev_line_list[args.column])]
    for line in in_fd:
        line_list = line.strip().split("\t")
        if line_list[0:3] == prev_line_list[0:3]:
            prev_line_list[-1] += float(line_list[args.column])
        else:
            prev_line_list[-1] = str(prev_line_list[-1])
            out_fd.write("\t".join(prev_line_list) + "\n")
            prev_line_list = line_list[0:3] + [float(line_list[args.column])]
    prev_line_list[-1] = str(prev_line_list[-1])
    out_fd.write("\t".join(prev_line_list) + "\n")

