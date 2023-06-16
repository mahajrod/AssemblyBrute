#!/usr/bin/env python
__author__ = "mahajrod"
"""
This file contains functions necessary for Snakemake file
"""
import os
import yaml
import shutil
from copy import deepcopy
from collections import OrderedDict
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path, PosixPath


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
    return  sorted(list(fastq_dir_path.glob("*{0}".format(fastq_extension))))


def copy_absent_entries(input_dictionary, output_dictionary):
    for entry in input_dictionary:
        if entry not in output_dictionary:
            output_dictionary[entry] = deepcopy(input_dictionary[entry])
        else:
            if not isinstance(output_dictionary[entry], Mapping): # check if existing entry is not dictionary or dictionary like
                continue # exit from recursion
            copy_absent_entries(input_dictionary[entry], output_dictionary[entry])


def detect_phasing_parameters(current_stage_parameters, phasing_stage, stage_separator=".."):
    parameter_list = current_stage_parameters.split(stage_separator)
    phasing_stage_coretools = []
    for settings in config["stage_coretools"][phasing_stage]:
        phasing_stage_coretools += config["stage_coretools"][phasing_stage][settings]
    stage_subparameters = None
    for entry in parameter_list:
        for tool in phasing_stage_coretools:
            if entry[:len(tool)] == tool:
                stage_subparameters = entry

    if stage_subparameters is None:
        raise ValueError("Impossible to detect phasing stage parameters for {0} and phasing stage {1}".format(current_stage_parameters,
                                                                                                              phasing_stage))
    for stage_parameters in stage_dict[phasing_stage]["parameters"].keys():
        if stage_subparameters in stage_parameters:
            return stage_parameters
    else:
        raise ValueError("Impossible to detect phasing stage parameters for {0} and phasing stage {1}".format(current_stage_parameters,
                                                                                                              phasing_stage))


def get_parameters_for_all_stages_in_chain(current_stage_parameters, stage_separator=".."):
    sub_parameter_list = current_stage_parameters.split(stage_separator)
    number_of_stages_in_chain = len(sub_parameter_list)
    chain_stage_dict = {}
    for index in range(0, number_of_stages_in_chain):
        parameters = stage_separator.join(sub_parameter_list[:index+1])
        for st in stage_dict:
            if "parameters" in stage_dict[st]:
                if parameters in stage_dict[st]["parameters"]:
                    stage = st
                    break
        else:
            raise ValueError("Impossible to recognize stage for parameters {0}".format(parameters))
        chain_stage_dict[stage] = parameters

    return chain_stage_dict


def get_input_assemblies(input_folder_path, ploidy, fasta_extention):
    print(input_folder_path, fasta_extention)
    fasta_filelist = list(input_folder_path.glob("*{0}".format(fasta_extention)))
    print(fasta_filelist)
    if len(fasta_filelist) != ploidy:
        raise ValueError("ERROR!!! Number of input fasta files ({0}) differs from ploidy ({1})!".format(len(fasta_filelist),
                                                                                                       ploidy))
    if ploidy == 1:
        return {"hap0": fasta_filelist[0]}
    else:
        fasta_dict = {}
        for hap in range(1, ploidy+1):
            haplotype = "hap{0}".format(hap)
            suffix = "hap{0}{1}".format(hap, fasta_extention)
            for filename in fasta_filelist:
                print(filename.name)
                if filename.name[:-len(suffix)] == suffix:

                    fasta_dict[haplotype] = filename.name
                    break
            else:
                raise ValueError("ERROR!!! Fasta file for haplotype hap{0} was not found!".format(hap))
        return fasta_dict
