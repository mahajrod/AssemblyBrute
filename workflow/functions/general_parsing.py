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
    for entry in parameter_list:
        for tool in phasing_stage_coretools:
            if entry[:len(tool)] == tool:
                print(entry)
                raise AssertionError("")
                return entry
    else:
        raise ValueError("Impossible to detect phasing stage parameters for {0} and phasing stage {1}".format(current_stage_parameters,
                                                                                                              phasing_stage))
