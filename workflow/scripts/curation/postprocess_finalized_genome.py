#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import glob
import sys
import shutil
import argparse

from pathlib import Path
from datetime import datetime
from functools import partial

from collections import OrderedDict

import pandas as pd

from RouToolPa.Parsers.Sequence import CollectionSequence
from RouToolPa.GeneralRoutines import FileRoutines
from ete3.tools.ete_build_lib.ordereddict import OrderedDict

DATE_FORMATS = ["%Y-%m-%d",
                "%d/%m/%Y",]

def parse_date(date_str):
    for fmt in DATE_FORMATS:
        try:
            return datetime.strptime(date_str, fmt)
        except ValueError:
            continue
    raise ValueError(f"ERROR!!! Incorrect date format: {date_str}. Use '%Y-%m-%d' or '%d/%m/%Y'")

def get_haplotype_files(prefix, ploidy, extension_list=[".fasta", ".fasta.gz", ".fa", ".fa.gz", ".fna", ".fna.gz"]):
    haplotype_suffix_list = ["hap0"] if ploidy == 1  else [f"hap{haplotype}" for haplotype in range(1, ploidy + 1)]
    haplotype_filename_dict = OrderedDict()
    for haplotype_suffix in haplotype_suffix_list:
        haplotype_filename_dict[haplotype_suffix] = []
        for extension in extension_list:
            haplotype_filename_dict[haplotype_suffix] += list(glob.glob(f"{prefix}*{haplotype_suffix}*{extension}"))

        if len(haplotype_filename_dict[haplotype_suffix]) > 1:
            raise ValueError(f"ERROR!!! Found more than one file with extensions: {','.join(extension_list)} "
                             f"for prefix = {prefix} and haplotype suffix = {haplotype_suffix}")
        elif len(haplotype_filename_dict[haplotype_suffix]) == 0:
            ValueError(f"ERROR!!! Found no files with extensions: {','.join(extension_list)} "
                       f"for prefix = {prefix} and haplotype suffix = {haplotype_suffix}")
        else:
            pass

    return haplotype_filename_dict

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
            description +=" [topology=circular]"
        if circular_mtdna:
            description+= " [completeness=complete]"

    return pd.NA if description == "" else description

def get_chromosome_name(scaffold_id):
    if ("aut" in scaffold_id) or ("chr" in scaffold_id): # chromosomes/autosomes
        return  scaffold_id.split("_")[0][3:]
    elif "cand" in scaffold_id: # candidate chromosomes, usually used for candidate sex chromosomes
        return  scaffold_id.split("_")[0]
    else:
        return None

def check_if_main_chr_scaffold(scaffold_id):
    if ("aut" in scaffold_id) or ("chr" in scaffold_id) or ("cand" in scaffold_id):
        if "unloc" in scaffold_id:
            return "no"
        return "yes"
    else:
        return "no"

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_prefix", action="store", dest="input_prefix", required=True,
                    help="Prefix of the haplotype fasta files. "
                         "Input files must follow <prefix>*hap<N>*<extension> or <prefix>*hap<N>*<extension> templates. "
                         "N=0 for haploid assemblies, and N=1..n for assembly with n > 1 haplotypes."
                         "Extention might be any of '.fasta', '.fa', '.fna'. Might be gzipped."
                         "* in the templates stands for any number of any symbols. Required")
parser.add_argument("-d", "--output_dir", action="store", dest="output_dir", default="./",
                    help="Directory to write output files. Default: current directory")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of the output files.Usage of versioned ToLID ids is highly recommended. "
                         "For example, bOtiPha1.1 for the first version of the genome. "
                         "Register your individual at https://id.tol.sanger.ac.uk/. Required.")
parser.add_argument("-l", "--ploidy", action="store", dest="ploidy", default=2, type=int,
                    help="Ploidy of the input assembly. Default: 2")
parser.add_argument("-n", "--naming_style", action="store", dest="naming_style", default="VGP",
                    help="Naming style of the output files. Allowed: 'VGP' (i.e. bOtiNob1.1.hap1.cur20260331.fasta.gz). Default: 'VGP'")
parser.add_argument("-e", "--date", action="store", dest="date",
                    help="Date to use for VGP style naming. Use 'YYYY-MM-DD' or 'DD/MM/YYYY'. If not set, the current date will be used")
parser.add_argument("-k", "--keep", action="store_true", dest="keep", default=False,
                    help="Keep intermediate files in case of exception. Default: False")
parser.add_argument("-g", "--gzip", action="store_true", dest="gzip", default=False,
                    help="Gzip output. Default: False")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Threads to use to gzip output. Default: 1")
parser.add_argument("-a","--not_complete_mtdna", action="store_true", dest="not_complete_mtdna",
                    default=False,
                    help="Mitochondrial DNA (scaffold is detected by 'mtDNA' name) is not complete. "
                         "Default: False, i.e, by default the script assumes a complete mtDNA.")
parser.add_argument("-b","--not_circular_mtdna", action="store_true", dest="not_circular_mtdna",
                    default=False,
                    help="Mitochondrial DNA (scaffold is detected by 'mtDNA' name) was not circularized during assembly. "
                         "Default: False, i.e. , by default the script assumes a circularized mtDNA ")

args = parser.parse_args()

curation_date = parse_date(args.date) if args.date else datetime.today()

# Recognize input fasta files
hap_filename_dict = get_haplotype_files(args.input_prefix, args.ploidy)
# get output directory
out_dir_path = Path(args.output_dir)

if args.naming_style == "VGP":
    current_date_str = curation_date.strftime('%Y%m%d')
    vgp_output_suffix = f"cur{current_date_str}"
else:
    raise ValueError("ERROR!!! Unrecognized naming style of the output! Check help for --naming_style option.")

tmp_dir_path = out_dir_path/ f"tmp_{os.urandom(30).hex()}"
os.mkdir(tmp_dir_path)

try:
    output_filepath_dict = OrderedDict()
    tmp_filepath_dict = OrderedDict()

    for haplotype in hap_filename_dict:
        if args.naming_style == "VGP":
            main_output_prefix = f"{args.output_prefix}.{haplotype}.{vgp_output_suffix}"

        tmp_filename_dict = {}
        output_filename_dict = {}
        for ext in "len", "description", "reorderlist", "fasta", "whitelist", "orderlist", "chromosomes":
            tmp_filename_dict[ext] =  tmp_dir_path / f"{main_output_prefix}.{ext}"
            output_filename_dict[ext] = out_dir_path / f"{main_output_prefix}.{ext}"

        # read fasta
        fasta_collection = CollectionSequence(in_file=hap_filename_dict[haplotype], get_stats=True)
        fasta_collection.seq_lengths.to_csv(tmp_filename_dict["len"], header=False, sep="\t")

        # create NCBI description
        len_df = pd.read_csv(tmp_filename_dict["len"], sep='\t', header=None, names=["scaffold_id", "length"], )
        len_df["description"] = len_df["scaffold_id"].apply(partial(create_description,
                                                                    complete_mtdna=not args.not_complete_mtdna,
                                                                    circular_mtdna=not args.not_circular_mtdna))

        len_df[["scaffold_id", "description"]][~len_df["description"].isna()].to_csv(tmp_filename_dict["description"],
                                                                                     sep="\t", header=False,
                                                                                     index=False)
        fasta_collection.description = len_df[["scaffold_id", "description"]].set_index("scaffold_id")["description"].to_dict()

        # create reorder list
        os.system(f" cut -f 1 {str(tmp_filename_dict['len'])} | grep -P 'aut|chr|cand|mtDNA' "
                  f"  |  sort -V > {str(tmp_filename_dict['reorderlist'])}") # replace by python code later, or external python  library - natsort
        os.system(f" grep -v -P 'mtDNA|unloc' {str(tmp_filename_dict['reorderlist'])} > {str(tmp_filename_dict['whitelist'])}")
        os.system(f" cp {str(tmp_filename_dict['whitelist'])}  {str(tmp_filename_dict['orderlist'])}")

        reorder_df = pd.read_csv(tmp_filename_dict['reorderlist'], sep="\t", header=None, names=["scaffold"])

        # reorder fasta
        fasta_collection.reorder_records(by="orderlist", orderlist=reorder_df['scaffold'], in_place=True)

        # create Sanger-style chromosome file
        reorder_df = reorder_df[reorder_df['scaffold'] != 'mtDNA']
        reorder_df['chromosome'] = reorder_df['scaffold'].apply(get_chromosome_name)
        reorder_df['main_chromosome_part'] = reorder_df['scaffold'].apply(check_if_main_chr_scaffold)
        reorder_df.to_csv(tmp_filename_dict["chromosomes"], sep=",", header=False,index=False)

        # write fasta
        fasta_collection.write(str(tmp_filename_dict['fasta']))

        for ext in tmp_filename_dict:
            shutil.move(tmp_filename_dict[ext], output_filename_dict[ext])

        if args.gzip:
            os.system(f"pigz -p {args.threads} {str(output_filename_dict['fasta'])}")

    if not args.keep:
        shutil.rmtree(tmp_dir_path)

except Exception as e:
    # remove temporary files in case of fail
    if not args.keep:
        shutil.rmtree(tmp_dir_path)
    raise

