#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
import numpy as np
import pysam

from RouToolPa.GeneralRoutines import FileRoutines

parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    type=lambda s: FileRoutines.metaopen(s, "r"),
                    help="Input bam file. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    type=lambda s: FileRoutines.metaopen(s, "w"),
                    help="Prefix of output files")
parser.add_argument("-l", "--length", action="store", dest="length", default=200, type=int,
                    help="Length of regions at the ends to assess")
parser.add_argument("-q", "--max_quality", action="store", dest="max_quality", default=100, type=int,
                    help="Maximal base quality. If exceeded base quality will be reduced to this value")
args = parser.parse_args()


bam = pysam.AlignmentFile(args.input, "rb")

start_nucleotide_df = pd.DataFrame(0, columns=["A", "C", "G", "T", "N"],
                                   index=range(0, args.length), dtype='Int64')
end_nucleotide_df = pd.DataFrame(0, columns=["A", "C", "G", "T", "N"],
                                 index=range(-1, -args.length-1, -1), dtype='Int64')

start_qc_df = pd.DataFrame(0, columns=range(0, args.max_quality+1),
                           index=range(0, args.length), dtype='Int64')
end_qc_df = pd.DataFrame(0, columns=range(0, args.max_quality+1),
                         index=range(-1, -args.length-1, -1), dtype='Int64')
print(end_qc_df)
min_read_length = 2 * args.length
read_number = 0
long_read_number = 0
for read in bam.fetch(until_eof=True):
    read_number += 1
    if read.query_length > min_read_length:
        long_read_number += 1
        for base_index in range(0, args.length):
            start_nucleotide_df.loc[base_index, read.query_sequence[base_index]] += 1
            start_qc_df.loc[base_index, min([read.query_qualities[base_index], args.max_quality])] += 1

        for base_index in range(-1, -args.length - 1, -1):
            end_nucleotide_df.loc[base_index, read.query_sequence[base_index]] += 1
            end_qc_df.loc[base_index, min([read.query_qualities[base_index], args.max_quality])] += 1

    try:
        print(read.query_name,
              read.query_length,
              np.mean(read.query_qualities),
              np.median(read.query_qualities),
              type(read.get_tags()))
    except:
        print(read)
