#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
from pathlib import Path


def convert_interval_format_to_bed(input, output):
    with open(input, "r") as in_fd:
        with open(output, "w") as out_fd:
            for line in in_fd:
                if line[0] == ">":
                    seq_id = line.strip().split()[0][1:]
                else:
                    start, stop = line.strip().split(" - ")
                    bed_string = "{0}\t{1}\t{2}\n".format(seq_id, start, stop)
                    out_fd.write(bed_string)


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input file with repeats in format produced by Windowmasker. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output",
                    default=sys.stdout,
                    help="Output bedfile file - default: stdout")

args = parser.parse_args()

convert_interval_format_to_bed(args.input, args.output)
