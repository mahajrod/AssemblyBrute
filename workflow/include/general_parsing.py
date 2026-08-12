#!/usr/bin/env python
__author__ = "mahajrod"
"""
This file contains functions necessary for Snakemake file
"""
from copy import deepcopy
from pathlib import Path, PosixPath


def detect_phasing_parameters(current_stage_parameters, phasing_stage, stage_separator=".."):
    parameter_list = current_stage_parameters.split(stage_separator)

    stage_subparameters = None
    stage_subparameters_index = None
    for index in range(0, len(parameter_list)):
        for tool in config["stage_coretools"][phasing_stage]:
            if parameter_list[index][:len(tool)] == tool:
                stage_subparameters = parameter_list[index]
                stage_subparameters_index = index

    if stage_subparameters is None:
        raise ValueError("Impossible to detect phasing stage parameters for {0} and phasing stage {1}".format(current_stage_parameters,
                                                                                                              phasing_stage))
    return stage_separator.join(parameter_list[:stage_subparameters_index+1])


def get_parameters_for_all_stages_in_chain(current_stage_parameters, stage_separator=".."):
    sub_parameter_list = current_stage_parameters.split(stage_separator)
    number_of_stages_in_chain = len(sub_parameter_list)
    chain_stage_dict = {}
    for index in range(0, number_of_stages_in_chain):
        parameters = stage_separator.join(sub_parameter_list[:index+1])
        for st in stage_dict:
            if hasattr(stage_dict[st], "parameters"):
                if parameters in stage_dict[st].parameters:
                    stage = st
                    break
        else:
            raise ValueError("Impossible to recognize stage for parameters {0}".format(parameters))
        chain_stage_dict[stage] = parameters

    return chain_stage_dict


def get_number_of_stages_in_chain(wildcards):
    return len(get_parameters_for_all_stages_in_chain(wildcards.parameters))

def detect_input_type(datatype, datatype_dir):
    datatype_dir_path = datatype_dir if isinstance(datatype_dir, PosixPath) else Path(datatype_dir)
    input_filedict = {}
    for allowed_input_type in config["data"][datatype]:
        filedict = {}
        input_dir_path = datatype_dir_path / allowed_input_type
        print("AAAAAAAAAAAA")
        print(datatype)
        print(allowed_input_type)
        print(config["data"][datatype][allowed_input_type])
        print(config["data"][datatype][allowed_input_type]["allowed_in_exts"])
        for extension in config["data"][datatype][allowed_input_type]["allowed_in_exts"]:
            files = sorted(list(input_dir_path.glob("*{0}".format(extension))))
            if files:
                filedict[extension] = deepcopy(files)
        if len(filedict) > 1:
            raise ValueError("ERROR!!! Input files for {0} data have different extensions: {1}. ".format(datatype, ",".join(filedict.keys())) +
                             "It might be a sign of incorrect data. Rename files to have same extension. " +
                             "Allowed extensions: {0}".format(",".join(config["data"][datatype][allowed_input_type]["allowed_in_exts"])))
        if filedict:
            input_filedict[allowed_input_type] = deepcopy(filedict)

    if len(input_filedict) > 1:
        raise ValueError("ERROR!!! Input files for {0} data are of different type: {1} ".format(datatype,
                                                                                                ",".join(input_filedict.keys())) +
                         "It might be a sign of incorrect data. Convert all data to the same format."
                         "Allowed formats: {0}".format(",".join(config["data"][datatype])))
    if len(input_filedict) == 0:
        raise ValueError(f"ERROR!!! Input files for {datatype} data are absent! Add data or remove it from the main config file...")

    return input_filedict

def parse_stage_parameters_from_path(path):
    string_list = str(path).split("/") # in case if string is nota str, but a Path

    candidate_stage_list = []
    candidate_stage_index_list = []
    for allowed_stage in config["allowed_stage_list"]:
        for index in range(0, len(string_list)):
            if string_list[index] == allowed_stage:
                candidate_stage_list.append(allowed_stage)
                candidate_stage_index_list.append(index)
    if len(candidate_stage_list) > 1:
        raise ValueError(f"ERROR!!! Impossible to detect a stage for path: {path}. Candidate stages: {', '.join(candidate_stage_list)} ")
    elif len(candidate_stage_list) == 0:
        raise ValueError(f"ERROR!!! Impossible to detect a stage for path: {path}. No allowed stage was found in it.")
    else:
        stage = candidate_stage_list[0]
        element_index = candidate_stage_index_list[0]

    parameters = string_list[element_index + 1]

    prev_parameters = None if ".." not in parameters else "..".join(parameters.split("..")[:-1])

    result_dict = {"stage": stage,
                   "parameters": parameters,
                   "prev_parameters": prev_parameters}

    return result_dict

def get_prev_stage_parameters(parameters):
    return "..".join(parameters.split("..")[:-1])

def get_relative_path(target: str or Path, link: str or Path):
    link_path = Path(link)

    return Path(os.path.relpath(target, start=link_path.parent))

def create_relative_link(target: str or Path, link: str or Path):
    link_path = Path(link)
    if link_path.exists() or link_path.is_symlink():
        link_path.unlink()
    relative_path = Path(os.path.relpath(target, start=link_path.parent))
    link_path.symlink_to(relative_path)
    return relative_path