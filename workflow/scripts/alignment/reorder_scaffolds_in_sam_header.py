#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import bz2
import gzip
import argparse
import pandas as pd

if sys.version_info[0] == 3:
    from io import TextIOWrapper as file

def metaopen(filename, flags, buffering=None, compresslevel=5):
    if not isinstance(filename, str): # or isinstance(filename, gzip.GzipFile) or isinstance(filename, bz2.BZ2File):
        if isinstance(filename, file):
            return filename
        else:
            raise ValueError("ERROR!!! Not file object or str: {}".format(str(filename)))
    elif filename[-3:] == ".gz":
        return gzip.open(filename, flags + ("t" if "b" not in flags else ""), compresslevel=compresslevel)
    elif filename[-4:] == ".bz2":
        return bz2.open(filename, flags + ("t" if "b" not in flags else ""), compresslevel=compresslevel)
    else:
        if buffering is not None:
            return open(filename, flags, buffering=buffering)
        else:
            return open(filename, flags)

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input file with sam file. It should not contain the records. Default: stdin")
parser.add_argument("-r", "--orderlist", action="store", dest="orderlist", required=True,
                    type=lambda s: pd.read_csv(s, sep="\t", header=None).squeeze(),
                    help="Input orderlist file. Required.")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file. Default: stdout")

args = parser.parse_args()

prefix_lines_list = []
scaffold_line_df = []
suffix_lines_list = []

with metaopen(args.input, "r") as in_fd:
    for line in in_fd:
        if line[:3] != "@SQ":
            prefix_lines_list.append(line)
        else:
            scaffold_line_df.append(line.strip().split("\t"))
            break
    for line in in_fd:
        if line[:3] == "@SQ": # all SQ lines should form a single block, so it is safe to read SQ lines and lines with different tag following SQ block
            scaffold_line_df.append(line.strip().split("\t"))
        else:
            suffix_lines_list.append(line)

scaffold_line_df  = pd.DataFrame.from_records(scaffold_line_df, columns=["tag","scaf_tag","len_tag"])
scaffold_line_df["scaffold"] = scaffold_line_df["scaf_tag"].apply(lambda s: s[3:])
scaffold_line_df = scaffold_line_df.set_index("scaffold")

reordered_df = scaffold_line_df.loc[args.orderlist]
# note that agp_1 file reported by PretextView doesn't contain scaffolds smaller than pixel length.
# So if you use orderlist extracted from such agp_1 file, remaining_df has high chance to be nonempty.
# it is fixed by running CorrectAGP script, but for separation of curation units it is not necessary.
remaining_df = scaffold_line_df.loc[~scaffold_line_df.index.isin(args.orderlist)]
if not remaining_df.empty:
    reordered_df = pd.concat([reordered_df, remaining_df])

with metaopen(args.output, "w") as out_fd:
    for line in prefix_lines_list:
        out_fd.write(line)
    reordered_df.to_csv(out_fd, sep="\t", header=False, index=False)
    for line in suffix_lines_list:
        out_fd.write(line)
