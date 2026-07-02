#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd

from RouToolPa.Parsers.Sequence import CollectionSequence

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--fasta", action="store", dest="fasta", default=sys.stdin,
                    help="Input fasta. Default: stdin")
parser.add_argument("-d", "--description_file", action="store", dest="description",
                    required=True,
                    help="Two-column (scaffold_id and description) file with description. Required")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output fasta. Default: stdout")

args = parser.parse_args()

description_dict = pd.read_csv(args.description, sep="\t", header=None, index_col=0,
                               names=["scaffold_id", "description"])["description"].to_dict()
#print(description_dict)
sequence_collection = CollectionSequence(in_file=args.fasta)

sequence_collection.description = description_dict
sequence_collection.write(args.output)
