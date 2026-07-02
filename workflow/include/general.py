#!/usr/bin/env python
__author__ = "mahajrod"
"""
This file contains functions necessary for Snakemake file
"""

from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path, PosixPath
from collections.abc import Iterable

def p_distance(seq_a, seq_b, seq_len):
    dist = 0
    for i in range(0, seq_len):
        if seq_a[i] != seq_b[i]:
            dist += 1
    return dist


def get_common_prefix_ans_suffixes(seq_a, seq_b):
    seq_a_len = len(seq_a)
    seq_b_len = len(seq_b)

    prefix = ""
    for i in range(0, min(seq_a_len, seq_b_len)):
        if seq_a[i] != seq_b[i]:
           return prefix, seq_a[i:], seq_b[i:]
        prefix += seq_a[i]
    return prefix, "", ""


def convert_posixpath2str_in_dict(dictionary):
    output_dictionary = deepcopy(dictionary)
    for entry in output_dictionary:
        if isinstance(output_dictionary[entry], PosixPath):
            output_dictionary[entry] = str(output_dictionary[entry])
        else:
            if not isinstance(output_dictionary[entry], Mapping): # check if existing entry is not dictionary or dictionary like
                continue # exit from recursion
            output_dictionary[entry] = convert_posixpath2str_in_dict(output_dictionary[entry])

    return output_dictionary


def find_cmap(bionano_dir, cmap_extension=".cmap"): # TODO: modify when input for bionano will be clear
    bionano_dir_path = bionano_dir if isinstance(bionano_dir, PosixPath) else Path(bionano_dir)
    cmap_list = list(bionano_dir_path.glob("*{0}".format(cmap_extension)))
    if len(cmap_list) > 1:
        raise ValueError("ERROR!!! More than one cmap file was found: {0}".format(", ".join(list(map(str, cmap_list)))
                                                                                  )
                         )
    return cmap_list[0]


def find_fastqs(fastq_dir, fastq_extension=".fastq.gz"):
    fastq_dir_path = fastq_dir if isinstance(fastq_dir, PosixPath) else Path(fastq_dir)
    return sorted(list(fastq_dir_path.glob("*{0}".format(fastq_extension))))


def find_fastas(fasta_dir, fasta_extension=".fasta.gz"):
    fasta_dir_path = fasta_dir if isinstance(fasta_dir, PosixPath) else Path(fasta_dir)
    return sorted(list(fasta_dir_path.glob("*{0}".format(fasta_extension))))


def find_bams(bam_dir, bam_extension=".bam"):
    bam_dir_path = bam_dir if isinstance(bam_dir, PosixPath) else Path(bam_dir)
    return sorted(list(bam_dir_path.glob("*{0}".format(bam_extension))))


def copy_absent_entries(input_dictionary, output_dictionary):
    for entry in input_dictionary:
        if entry not in output_dictionary:
            output_dictionary[entry] = deepcopy(input_dictionary[entry])
        else:
            if not isinstance(output_dictionary[entry], Mapping): # check if existing entry is not dictionary or dictionary like
                continue # exit from recursion
            copy_absent_entries(input_dictionary[entry], output_dictionary[entry])

def log_formatted_nested_dict(dictionary, logger, indent_type, initial_indent_level, prefix=""):

    for subdict in dictionary:
        if isinstance(dictionary[subdict], str) or not(isinstance(dictionary[subdict], Iterable)):
            logger(prefix + indent_type * initial_indent_level + "- " + str(subdict) + ":" + indent_type + str(dictionary[subdict]))
            continue
        elif isinstance(dictionary[subdict], (list, set)):
            logger(prefix + indent_type * initial_indent_level + "- " + str(subdict) + ":")
            for element in dictionary[subdict]:
                logger(prefix + indent_type * (initial_indent_level + 1) + "- " + str(element))
            continue
        logger(prefix + indent_type * initial_indent_level + "- " + str(subdict) + ":")
        log_formatted_nested_dict(dictionary[subdict], logger, indent_type, initial_indent_level+1, prefix=prefix)


def get_input_assemblies(input_folder_path, ploidy, fasta_extention):
    fasta_filelist = list(input_folder_path.glob("*{0}".format(fasta_extention)))
    if len(fasta_filelist) != ploidy:
        raise ValueError("ERROR!!! Number of input fasta files ({0}) differs from ploidy ({1})!".format(len(fasta_filelist),
                                                                                                       ploidy))
    if ploidy == 1:
        return {"hap0": fasta_filelist[0].name}
    else:
        fasta_dict = {}
        for hap in range(1, ploidy+1):
            haplotype = "hap{0}".format(hap)
            suffix = "hap{0}{1}".format(hap, fasta_extention)
            for filename in fasta_filelist:
                if filename.name[-len(suffix):] == suffix:
                    fasta_dict[haplotype] = filename.name
                    break
            else:
                raise ValueError("ERROR!!! Fasta file for haplotype hap{0} was not found!".format(hap))
        return fasta_dict

