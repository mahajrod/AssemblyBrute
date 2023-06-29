#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
from RouToolPa.Parsers.Sequence import CollectionSequence
from RouToolPa.GeneralRoutines import FileRoutines

parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input fasta file. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output bedgraph file. Default: stdout")
parser.add_argument("-w", "--window_size", action="store", dest="window_size", default=100000, type=int,
                    help="Size of the windows Default: 100000")
parser.add_argument("-s", "--window_step", action="store", dest="window_step", default=None, type=int,
                    help="Step of the sliding windows. Default: window size, i.e windows are staking")

args = parser.parse_args()

if args.window_step is None:
    args.window_step = args.window_size

fasta_collection = CollectionSequence(in_file=args.input)

with FileRoutines.metaopen(args.output, "w") as out_fd:
    for scaffold in fasta_collection.records:
        sequence_length = len(fasta_collection.records[scaffold])
        if sequence_length < args.window_size:
            continue
        number_of_steps = ((sequence_length - args.window_size) // args.window_step) + 1

        for window_id in range(0, number_of_steps):
            start = window_id * args.window_step
            end = window_id * args.window_step + args.window_size
            out_fd.write("{0}\t{1}\t{2}\t{3}\n".format(scaffold, start, end,
                                                       (fasta_collection.records[scaffold][start:end].count("G") +
                                                        fasta_collection.records[scaffold][start:end].count("C") +
                                                        fasta_collection.records[scaffold][start:end].count("g") +
                                                        fasta_collection.records[scaffold][start:end].count("c")) / args.window_size))




