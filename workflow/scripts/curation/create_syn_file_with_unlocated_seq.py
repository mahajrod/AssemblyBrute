#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import signal
import argparse

from functools import partial
import pandas as pd

from RouToolPa.Collections.General import SynDict, IdList

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input file with ids to create syns for. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output TAB-separated file with new synonyms. Column 0 contains original ids, column 1 - new synonyms, respectively. Default: stdout")
parser.add_argument("-c", "--chr_syn_file", action="store", dest="chr_syn_file", required=True,
                    help="File with chromosome syns. Required")
parser.add_argument("-s", "--id_separator", action="store", dest="id_separator", default="_",
                    help="Separator used to split chr and unloc id. Default: '_'")
parser.add_argument("-u", "--unloc_label", action="store", dest="unloc_label", default="unloc",
                    help="Label of unlocalized sequences. Default: 'unloc'")
parser.add_argument("-k", "--syn_file_key_column", action="store", dest="syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in chromosome synonym file. Default: 0")
parser.add_argument("-v", "--syn_file_value_column", action="store", dest="syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in chromosome synonym file. Default: 1")

args = parser.parse_args()

chr_syn_dict = SynDict(filename=args.chr_syn_file,
                       key_index=args.syn_file_key_column,
                       value_index=args.syn_file_value_column)

combined_separator = args.id_separator + args.unloc_label

def split_orig_id(orig_id, combo_separator, unloc_label):
    tmp = orig_id.split(combo_separator)
    if len(tmp) == 1:
        return pd.Series([tmp[0], pd.NA])
    elif len(tmp) == 2:
        return pd.Series([tmp[0], unloc_label + tmp[1]])
    else:
        raise ValueError(f"ERROR!!! Combined separator '{combined_separator}' is preszent more than one time in the original scaffold id '{orig_id}'")

def create_renamed_id(row, id_separator):
    if pd.isna(row.iloc[1]):
        return row.iloc[0]
    else:
        return row.iloc[0] + id_separator + row.iloc[1]

id_df = pd.read_csv(args.input, header=None, dtype=str, names=["orig_ids"])
split_ids_df = id_df["orig_ids"].apply(partial(split_orig_id, combo_separator=combined_separator, unloc_label=args.unloc_label))
split_ids_df.columns = pd.Index(["chr_id", "unloc_id"])
id_df = pd.concat([id_df, split_ids_df], axis=1).set_index("orig_ids")
id_df["renamed_chr_id"] = id_df["chr_id"].replace(chr_syn_dict)
id_df["renamed_id"] = id_df[["renamed_chr_id", "unloc_id"]].apply(partial(create_renamed_id, id_separator=args.id_separator), axis=1)
id_df[["renamed_id"]].to_csv(args.output, sep="\t", header=False, index=True)
