#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
import pandas as pd

from functools import partial

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input yggdrasil len file. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file with description column for submission. Default: stdout")
parser.add_argument("-a","--not_complete_mtdna", action="store_true", dest="not_complete_mtdna",
                    default=False,
                    help="Mitochondrial DNA (scaffold is detected by 'mtDNA' name) is not complete. "
                         "Default: False, i.e, by default the script assumes a complete mtDNA.")
parser.add_argument("-b","--not_circular_mtdna", action="store_true", dest="not_circular_mtdna",
                    default=False,
                    help="Mitochondrial DNA (scaffold is detected by 'mtDNA' name) was not circularized during assembly. "
                         "Default: False, i.e. , by default the script assumes a circularized mtDNA ")

args = parser.parse_args()

def create_description(scaffold_id, complete_mtdna=True, circular_mtdna=True):
    description = ""
    if ("aut" in scaffold_id) or ("chr" in scaffold_id): # chromosomes/autosomes
        chr_id = scaffold_id.split("_")[0][3:]
        description = "[chromosome={0}]".format(chr_id)
        if "unloc" not in scaffold_id:
            description += " [location=chromosome]"
    elif "cand" in scaffold_id: # candidate chromosomes, usually used for candidate sex chromosomes
        chr_id = scaffold_id.split("_")[0]
        description = "[chromosome={0}]".format(chr_id)
        if "unloc" not in scaffold_id:
            description += " [location=chromosome]"
    elif scaffold_id == "mtDNA":
        description += "[location=mitochondrion]"
        if complete_mtdna:
            description += " [topology=circular]"
        if circular_mtdna:
            description += " [completeness=complete]"

    return pd.NA if description == "" else description

len_df = pd.read_csv(args.input, sep='\t', header=None, names=["scaffold_id", "length"],)
len_df["description"] = len_df["scaffold_id"].apply(partial(create_description,
                                                            complete_mtdna=not args.not_complete_mtdna,
                                                            circular_mtdna=not args.not_circular_mtdna))

len_df[["scaffold_id", "description"]][~len_df["description"].isna()].to_csv(args.output, sep="\t", header=False, index=False)
