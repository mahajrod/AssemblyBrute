#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd
from RouToolPa.Parsers.Sequence import CollectionSequence


parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input file with candidate microchromosome scaffolds. Default: stdin")
parser.add_argument("-a", "--assembly", action="store", dest="assembly_file", default=None,
                    help="Assembly file. Default: not set")
parser.add_argument("-f", "--fasta", action="store", dest="fasta", default=None,
                    help="Fasta file. Default: not set")
parser.add_argument("-l", "--len_file", action="store", dest="len_file", required=True,
                    help="File with scaffold lengths")
parser.add_argument("-m", "--max_length", action="store", dest="max_length", default=10000000, type=int,
                    help="Maximal length of the candidate microchromosome scaffold. Default: 10'000'000")
parser.add_argument("-s", "--scaffold_id_prefix", action="store", dest="scaffold_id_prefix", default="scaffold_",
                    help="Prefix of the scaffold ids. Required to reorder assembly file. "
                         "Default: 'scaffold_' (for yahs generated files). For 3ddna files use 'HiC_scaffold_'")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")

args = parser.parse_args()

len_df = pd.read_csv(args.len_file, sep="\t", header=None, names=["scaffold_id", "length"], index_col="scaffold_id")

candidate_scaffold_df = pd.read_csv(args.input, sep="\t", names=["scaffold_id", "protein_number"], header=None)

final_scaffold_df = candidate_scaffold_df[candidate_scaffold_df["scaffold_id"].isin(len_df[len_df["length"] <= args.max_length].index)].set_index("scaffold_id", inplace=False)

final_scaffold_df["length"] = len_df.loc[final_scaffold_df.index]["length"]

final_scaffold_df.to_csv(args.output_prefix + ".reordered.candidates.microchromosomes.filtered.tsv", index=True,
                         header=True, sep="\t")

if args.assembly_file:
    final_scaffold_index_list = sorted(list(map(lambda s: int(s[len(args.scaffold_id_prefix):]),
                                         list(final_scaffold_df.index))))
    with open(args.assembly_file, "r") as assembly_fd, \
         open(args.output_prefix + ".reordered.assembly", "w") as out_fd:
        end_list = []
        for line in assembly_fd:
            if line[0] == ">":
                out_fd.write(line)
            else:
                if 1 in final_scaffold_index_list:
                    out_fd.write(line)
                else:
                    end_list.append(line)
                index = 2
                for scaf_line in assembly_fd:
                    if index in final_scaffold_index_list:
                        out_fd.write(scaf_line)
                    else:
                        end_list.append(scaf_line)
                    index += 1
        for line_string in end_list:
            out_fd.write(line_string)
if args.fasta:
    print("Writing fasta...")
    collection_fasta = CollectionSequence(in_file=args.fasta)
    scaffold_series = pd.Series(collection_fasta.records.keys())
    non_microchromosome_scaffolds = scaffold_series[~scaffold_series.isin(final_scaffold_df.index)]
    collection_fasta.write(args.output_prefix + ".reordered.fasta",
                           whitelist=list(final_scaffold_df.index) + list(non_microchromosome_scaffolds))
#with open(args.output_prefix + "", "w") as out_fd:

