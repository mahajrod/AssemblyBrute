#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
from RouToolPa.Parsers.Sequence import CollectionSequence
from RouToolPa.GeneralRoutines import FileRoutines

parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input sam file without read qualities, for example generated from alignment of "
                         "error-corrected pacbio reads. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output sam file. Default: stdout")
parser.add_argument("-q", "--quality", action="store", dest="quality", default=30, type=int,
                    help="Quality value to use. Default: 30")
parser.add_argument("-e", "--encoding_shift", action="store", dest="encoding_shift", default=33, type=int,
                    help="Shift of quality encoding. Default: '33', i.e phred33")

args = parser.parse_args()


quality_symbol = chr(args.quality + args.encoding_shift)

quality_string_length = 10**7
quality_string = quality_symbol * quality_string_length

with FileRoutines.metaopen(args.input, "r") as in_fd, FileRoutines.metaopen(args.output, "r") as out_fd:
    for line in in_fd:
        if line[0] == "@":
            out_fd.write(line)
            continue
        line_list = line.split("\t")
        read_length = len(line_list[9])
        line_list[10] = quality_string[:read_length] if read_length < quality_string_length else quality_symbol * read_length
        out_fd.write("\t".join(line_list))




