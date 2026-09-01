import os
import sys
from functools import partial
import yaml
import shutil
from copy import deepcopy
from collections import OrderedDict
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path, PosixPath
from numbers import Number # Abstract class for numeric types
import pandas as pd
import numpy as np
from pyexpat import features

from workflow.include.general import *
from workflow.include.logging import config_logger, TAB

from contextlib import redirect_stdout, redirect_stderr

#---- Include sections for functions and config ----
# TODO: reorganize most of the included code as library(ies) and normal imports
include: "workflow/include/option_parsing.py"
include: "workflow/include/general_parsing.py"
include: "workflow/include/resources.py"
include: "workflow/include/read_configs.py"
include: "workflow/include/stage.py"
#----
#---- Initialize AssemblyBrute logger ----
brute_logger = config_logger(config)
#----

#---- Set conda environment for singularity ----
singularity_conda_env_title = config["singularity_load_env"] if config["singularity_load_mode"] == "cluster" else "singularity"
singularity_conda_env = config["conda"][singularity_conda_env_title]["name"] if config["use_existing_envs"] else f"../../../{config['conda'][singularity_conda_env_title]['yaml']}"
#-----

#---- Initialization of path variables from config file ----
#-------- Initialization of path variables for output --------
config["out_dir"] = Path(config["out_dir"])
config["out_dict"] = {}
for dir_name in "log", "tmp", "config":
    if not (config["out_dir"] / dir_name).exists(): # verification is needed as recreation of this dirs (especially 'log') triggers rerun of multiple rules
        os.makedirs(config["out_dir"] / dir_name)
#--------
#-------- Initialization of path variables for input --------
input_dir_path = Path(config["input_dir"])
brute_logger.info("AssemblyBrute logger initialized...")

config["input_datatypes"] = config["input_datatypes"].split(",")
config["data_feature_dict"] = {feature: set() for feature in ["paired", "fastq", "fasta", "fastqc",
                                                              "long_read", "nanopore", "pacbio",
                                                              "genome_size", "variant_call", "gap_fill",
                                                              "kraken", "filter", "phasing", "pretext_coverage_track",
                                                              "pretext_per_hap_track"]}

brute_logger.info(f"Checking input files...")

config["data"] = deepcopy(config["data_parameters"])

for datatype in list(config["data"].keys()):
    if datatype not in config["input_datatypes"]:
        # remove absent datatypes
        config["data"].pop(datatype)
        brute_logger.dbg_scr(TAB + f"Skipping input {datatype} files (absent in the main config)...")
        continue
    if datatype in ["reference", "draft"]:
        #this datatypes are parsed later
        continue
    # parsing input filenames
    brute_logger.info(TAB + f"Checking input {datatype} files...")
    datatype_dir = input_dir_path / datatype
    input = detect_input_type(datatype, input_dir_path / datatype)
    input_format = list(input.keys())[0]
    config["data"][datatype] = config["data"][datatype][input_format]
    config["data"][datatype].pop("allowed_in_exts")
    config["data"][datatype]["in_ext"] = list(input[input_format].keys())[0]
    config["data"][datatype]["in_dir"] = input_dir_path / datatype / input_format
    config["data"][datatype]["in_file_list"] = input[input_format][config["data"][datatype]["in_ext"]]
    config["data"][datatype]["num_files"] = len(config["data"][datatype]["in_file_list"])
    config["data"][datatype]["file_prefix_list"] = list(map(lambda s: str(s.name)[:-len(config["data"][datatype]["in_ext"])],
                                                        config["data"][datatype]["in_file_list"]))

    if config["data"][datatype]["paired"]:
        config["data"][datatype]["pair_prefix_list"] = []
        config["data"][datatype]["conv_file_prefix_list"] = []
        if (config["data"][datatype]["num_files"] % 2) != 0:
            raise ValueError(f"ERROR!!! {datatype} fastq files seems to be unpaired or misrecognized")
        for forward, reverse in zip(config["data"][datatype]["in_file_list"][::2],
                                    config["data"][datatype]["in_file_list"][1::2]):
            if p_distance(str(forward), str(reverse), len(str(forward))) > 1:
                raise ValueError(f"ERROR!!! {datatype} forward and reverse read files differ by more than one symbol:\n\t{forward}\n\t{reverse}")
        config["data"][datatype]["in_fwd_sfx"] = set()
        config["data"][datatype]["in_rev_sfx"] = set()
        # detect pairprefix, forward_and_reverse_suffixes for paired data
        for forward_prefix, reverse_prefix in zip(config["data"][datatype]["file_prefix_list"][::2],
                                                  config["data"][datatype]["file_prefix_list"][1::2]):
            common_prefix, forward_suffix, reverse_suffix = get_common_prefix_ans_suffixes(forward_prefix, reverse_prefix)
            config["data"][datatype]["pair_prefix_list"].append(common_prefix)
            config["data"][datatype]["in_fwd_sfx"].add(forward_suffix)
            config["data"][datatype]["in_rev_sfx"].add(reverse_suffix)
        if (len(config["data"][datatype]["in_fwd_sfx"]) > 1) or (len(config["data"][datatype]["in_rev_sfx"]) > 1):
            raise ValueError(f"ERROR!!! Multiple different suffixes in {datatype} filenames!")

        config["data"][datatype]["in_fwd_sfx"] = list(config["data"][datatype]["in_fwd_sfx"])[0]
        config["data"][datatype]["in_rev_sfx"] = list(config["data"][datatype]["in_rev_sfx"])[0]

        for pairprefix in config["data"][datatype]["pair_prefix_list"]:
            config["data"][datatype]["conv_file_prefix_list"].append(pairprefix + config["data"][datatype]["conv_fwd_sfx"])
            config["data"][datatype]["conv_file_prefix_list"].append(pairprefix + config["data"][datatype]["conv_rev_sfx"])
    else: # register prefixes of files for se data to simplify dealing with wildcards
        config["data"][datatype]["pair_prefix_list"] = config["data"][datatype]["file_prefix_list"]
        config["data"][datatype]["conv_file_prefix_list"] = config["data"][datatype]["file_prefix_list"]

    # check datatype specific filtering requests
    config["data"][datatype]["filter"] = False if datatype not in config["data_filtering"] else (config["data"][datatype]["filter"] & True)

    # create output dirnames
    config["data"][datatype]["raw_dir"] = config["out_dir"] / "data" / datatype / "raw"
    config["data"][datatype]["trimmed_dir"] = config["out_dir"] / "data" / datatype / "trimmed"
    config["data"][datatype]["filtered_dir"] = config["out_dir"] / "data" / datatype / "filtered"
    config["data"][datatype]["final_dir"] = config["out_dir"] / "data" / datatype / "final"

    if config["data"][datatype]["conv_fmt"] == "fastq":
        config["data_feature_dict"]["fastq"].add(datatype)
    if config["data"][datatype]["conv_fmt"] == "fasta":
        config["data_feature_dict"]["fasta"].add(datatype)
    for feature in ("paired", "fastqc", "long_read", "nanopore", "pacbio", "genome_size", "variant_call", "gap_fill",
                    "kraken", "filter", "phasing", "pretext_coverage_track", "pretext_per_hap_track"):
        if config["data"][datatype][feature]:
            config["data_feature_dict"][feature].add(datatype)

    brute_logger.info(TAB * 2 + f"Input extension: {config['data'][datatype]['in_ext']}")
    brute_logger.info(TAB * 2 + f"Input files: {config['data'][datatype]['num_files']}")
    brute_logger.info(TAB * 2 + "Files:")
    for filepath in config["data"][datatype]["in_file_list"]:
        brute_logger.info(TAB * 3 + str(filepath))

config["ext_data_feature_dict"] = {}

config["ext_data"] = {}
if "ext_data" in config["input_datatypes"]: # parse data that will be used to create additional tracks for pretextview. It is not used for assembly itself
    brute_logger.info(TAB + f"Checking external data...")

    for datatype in list(config["data_parameters"].keys()):

        datatype_dir = input_dir_path / "ext_data" / datatype
        track_dir_list = []
        for track_dir in datatype_dir.glob("*"):
            if track_dir.is_dir():
                track_dir_list.append(track_dir)
        if not track_dir_list:
            brute_logger.dbg_scr(TAB * 2 + f"Checking external {datatype} data...")
            brute_logger.dbg_scr(TAB * 3 + f"External {datatype} data was not found...")
            continue
        brute_logger.info(TAB * 2 + f"Checking external{datatype} data...")
        config["ext_data"][datatype] = {}
        config["ext_data_feature_dict"][datatype] = {feature: set() for feature in ["paired", "fastq", "fasta", "fastqc",
                                                                                      "long_read", "nanopore", "pacbio",
                                                                                      "genome_size", "variant_call", "gap_fill",
                                                                                      "kraken", "filter", "phasing", "pretext_coverage_track",
                                                                                      "pretext_per_hap_track"]}
        for track_dir in track_dir_list:
            track_name = track_dir.name
            brute_logger.info(TAB * 3 + f"Found external dataset {track_name}...")
            input = detect_input_type(datatype, track_dir)

            input_format = list(input.keys())[0]

            config["ext_data"][datatype][track_name] = deepcopy(config["data_parameters"][datatype][input_format])
            config["ext_data"][datatype][track_name].pop("allowed_in_exts")
            config["ext_data"][datatype][track_name]["in_ext"] = list(input[input_format].keys())[0]
            config["ext_data"][datatype][track_name]["in_dir"] = track_dir / input_format
            config["ext_data"][datatype][track_name]["in_file_list"] = input[input_format][config["ext_data"][datatype][track_name]["in_ext"]]
            config["ext_data"][datatype][track_name]["num_files"] = len(config["data"][datatype]["in_file_list"])
            config["ext_data"][datatype][track_name]["file_prefix_list"] = list(map(lambda s: str(s.name)[:-len(config["ext_data"][datatype][track_name]["in_ext"])],
                                                                config["ext_data"][datatype][track_name]["in_file_list"]))

            if config["ext_data"][datatype][track_name]["paired"]:
                config["ext_data"][datatype][track_name]["pair_prefix_list"] = []
                config["ext_data"][datatype][track_name]["conv_file_prefix_list"] = []
                if (config["ext_data"][datatype][track_name]["num_files"] % 2) != 0:
                    raise ValueError(f"ERROR!!! {datatype} fastq files seems to be unpaired or misrecognized")
                for forward, reverse in zip(config["ext_data"][datatype][track_name]["in_file_list"][::2],
                                            config["ext_data"][datatype][track_name]["in_file_list"][1::2]):
                    if p_distance(str(forward), str(reverse), len(str(forward))) > 1:
                        raise ValueError(f"ERROR!!! {datatype} forward and reverse read files differ by more than one symbol:\n\t{forward}\n\t{reverse}")
                config["ext_data"][datatype][track_name]["in_fwd_sfx"] = set()
                config["ext_data"][datatype][track_name]["in_rev_sfx"] = set()
                # detect pairprefix, forward_and_reverse_suffixes for paired data
                for forward_prefix, reverse_prefix in zip(config["ext_data"][datatype][track_name]["file_prefix_list"][::2],
                                                          config["ext_data"][datatype][track_name]["file_prefix_list"][1::2]):
                    common_prefix, forward_suffix, reverse_suffix = get_common_prefix_ans_suffixes(forward_prefix, reverse_prefix)
                    config["ext_data"][datatype][track_name]["pair_prefix_list"].append(common_prefix)
                    config["ext_data"][datatype][track_name]["in_fwd_sfx"].add(forward_suffix)
                    config["ext_data"][datatype][track_name]["in_rev_sfx"].add(reverse_suffix)
                if (len(config["ext_data"][datatype][track_name]["in_fwd_sfx"]) > 1) or (len(config["ext_data"][datatype][track_name]["in_rev_sfx"]) > 1):
                    raise ValueError(f"ERROR!!! Multiple different suffixes in {datatype} filenames!")

                config["ext_data"][datatype][track_name]["in_fwd_sfx"] = list(config["ext_data"][datatype][track_name]["in_fwd_sfx"])[0]
                config["ext_data"][datatype][track_name]["in_rev_sfx"] = list(config["ext_data"][datatype][track_name]["in_rev_sfx"])[0]

                for pairprefix in config["ext_data"][datatype][track_name]["pair_prefix_list"]:
                    config["ext_data"][datatype][track_name]["conv_file_prefix_list"].append(pairprefix + config["ext_data"][datatype][track_name]["conv_fwd_sfx"])
                    config["ext_data"][datatype][track_name]["conv_file_prefix_list"].append(pairprefix + config["ext_data"][datatype][track_name]["conv_rev_sfx"])
            else: # register prefixes of files for se data to simplify dealing with wildcards
                config["ext_data"][datatype][track_name]["pair_prefix_list"] = config["ext_data"][datatype][track_name]["file_prefix_list"]
                config["ext_data"][datatype][track_name]["conv_file_prefix_list"] = config["ext_data"][datatype][track_name]["file_prefix_list"]

            # check datatype specific filtering requests
            config["ext_data"][datatype][track_name]["filter"] = False if datatype not in config["data_filtering"] else (config["ext_data"][datatype][track_name]["filter"] & True)

            # create output dirnames
            config["ext_data"][datatype][track_name]["raw_dir"] = config["out_dir"] / "ext_data" / datatype / track_name / "raw"
            config["ext_data"][datatype][track_name]["trimmed_dir"] = config["out_dir"] / "ext_data" / datatype / track_name / "trimmed"
            config["ext_data"][datatype][track_name]["filtered_dir"] = config["out_dir"] / "ext_data" / datatype / track_name / "filtered"
            config["ext_data"][datatype][track_name]["final_dir"] = config["out_dir"] / "ext_data" / datatype / track_name / "final"


            if config["ext_data"][datatype][track_name]["conv_fmt"] == "fastq":
                config["ext_data_feature_dict"][datatype]["fastq"].add(track_name)
            if config["ext_data"][datatype][track_name]["conv_fmt"] == "fasta":
                config["ext_data_feature_dict"][datatype]["fasta"].add(track_name)
            for feature in ("paired", "fastqc", "long_read", "nanopore", "pacbio", "genome_size", "variant_call", "gap_fill",
                            "kraken", "filter", "phasing", "pretext_coverage_track", "pretext_per_hap_track"):
                if config["ext_data"][datatype][track_name][feature]:
                    config["ext_data_feature_dict"][datatype][feature].add(track_name)

            brute_logger.info(TAB * 4 + f"Input extension: {config['ext_data'][datatype][track_name]['in_ext']}")
            brute_logger.info(TAB * 4 + f"Input files: {config['ext_data'][datatype][track_name]['num_files']}")
            brute_logger.info(TAB * 4 + "Files:")
            for filepath in config["ext_data"][datatype][track_name]["in_file_list"]:
                brute_logger.info(TAB * 5 + str(filepath))


if "reference" in config["data"]:
    brute_logger.info(TAB + f"Checking input reference files...")
    config["data"]["reference"]["in_dir"] = input_dir_path / "reference"
    config["data"]["reference"]["ref_dict"] = {}
    for filename in config["data"]["reference"]["in_dir"].glob("*"):
        if filename.is_dir():
            config["data"]["reference"]["ref_dict"][filename.name] = {}
    for genome in config["data"]["reference"]["ref_dict"]:
        brute_logger.info(TAB * 2 + f"Checking reference {genome}...")
        for filetype in "fasta", "syn", "whitelist", "orderlist":
            config["data"]["reference"]["ref_dict"][genome][filetype] = list((config["data"]["reference"]["in_dir"] / genome).glob(f"*.{filetype}"))

            if len(config["data"]["reference"]["ref_dict"][genome][filetype]) > 1:
                raise ValueError(f"ERROR!!! There is more than one {filetype} file for reference {genome}")

            config["data"]["reference"]["ref_dict"][genome][filetype] = config["data"]["reference"]["ref_dict"][genome][filetype][0]
            brute_logger.info(TAB * 3 + f"Detected {filetype}:")
            brute_logger.info(TAB * 4 + f"{config['data']['reference']['ref_dict'][genome][filetype]}")

if "draft" in config["data"]:
    brute_logger.info(TAB + f"Checking input draft files...")
    config["data"]["draft"]["in_dir"] = input_dir_path / "draft/fasta"
    config["data"]["draft"]["haplotypes"] = get_input_assemblies(input_dir_path / "draft/fasta", config["ploidy"], config["data"]["draft"]["fasta"]["allowed_in_ext"])
    brute_logger.info(TAB * 2  + f"Detected haplotypes ({len(config['data']['draft']['haplotypes'])}):")
    for haplotype in config["data"]["draft"]["haplotypes"]:
        brute_logger.info(TAB * 3  + str(config["data"]["draft"]["haplotypes"][haplotype]))

#------------------------------------------------------------------------------------------

brute_logger.info("")
brute_logger.low_dbg_scr("Data dict after input file check:")
log_formatted_nested_dict(config["data"], brute_logger.low_dbg_scr, TAB, initial_indent_level=1, prefix="")
brute_logger.info("")
brute_logger.low_dbg_scr("Data feature dict after input file check:")
log_formatted_nested_dict(config["data_feature_dict"], brute_logger.low_dbg_scr, TAB, initial_indent_level=1, prefix="")

#---- Parse AGP_1 file with curation units ----
# It should be an AGP_1 file generated by PretextView
candidate_agp_dir_path = input_dir_path / "candidate_chr/"

candidate_agp_filename = list(candidate_agp_dir_path.glob("*.agp_1"))

if len(candidate_agp_filename) > 1:
    raise ValueError(f"ERROR!!! More than one agp file was detected in folder {str(candidate_agp_dir_path.name)}!")
elif len(candidate_agp_filename) == 1:
    candidate_agp_filename = candidate_agp_filename[0]

    brute_logger.info("")
    brute_logger.info("Detected AGP file with curation units:")
    brute_logger.info(TAB + str(candidate_agp_filename))

    candidate_output_dir = config["out_dir"] / "data/candidate_chr/"
    if not candidate_output_dir.exists():
        os.system(f" mkdir -p {str(candidate_output_dir)}")
    candidate_output_prefix = candidate_output_dir / "candidate"
    agp_df = pd.read_csv(candidate_agp_filename, sep="\t", header=None,
                            names=["scaffold_id", "start", "end", "part_number", "part_type",
                                   "part_id/gap_length", "part_start/gap_type",
                                   "part_end/linkage", "orientation/evidence", "comment"],
                            comment="#",
                            index_col="scaffold_id", usecols=[0,1,2,3,4,5,6,7,8,9])

    all_contig_series = agp_df[agp_df["part_type"] != "U"]["part_id/gap_length"]
    chr_component_series = agp_df[agp_df["comment"] == "Painted"]["part_id/gap_length"]

    chr_component_series.to_csv(f"{candidate_output_prefix}.all_chr.components.ids",sep="\t",header=False,index=False)
    candidate_chr_id_list = list(chr_component_series.index)
    brute_logger.info("Curation units:")
    for curation_unit in chr_component_series.index.unique():
        brute_logger.info(TAB + f"- {curation_unit}")
        for scaffold_id in chr_component_series[curation_unit]:
            brute_logger.dbg_scr(TAB * 2 + "- " + scaffold_id)
        chr_component_series[[curation_unit]].to_csv(f"{candidate_output_prefix}.{curation_unit}.components.ids",sep="\t",header=False,index=False)
        chr_black_list_series = chr_component_series[~chr_component_series.isin(chr_component_series[[curation_unit]])]
        chr_black_list_series.to_csv(f"{candidate_output_prefix}.{curation_unit}.pretext.blacklist",sep="\t",header=False,index=False)
else:
    candidate_chr_id_list = []
    candidate_agp_filename = []
#----

#---- Initialize tool parameters ----
#logging.info("Initializing tool parameters...")
#check if custom restriction sites were provided:

if config["custom_enzyme_set"] is not None:
    config["hic_enzyme_set"] = "custom"

    if config["custom_enzyme_set_is_no_motif"]: # register custom enzyme as producing no ligation motives
        config["no_motif_enzyme_sets"].append("custom")

if config["parameter_set"] not in config["parameters"]:
    raise ValueError(f"Error!!! Unknown set of tool parameters: {config['parameter_set']}")

copy_absent_entries(config["parameters"]["default"], config["parameters"][config["parameter_set"]]) # set default values for options absent in  "parameter_set"

for key in list(config["parameters"].keys()): # remove unused sets of parameters
    if key != config["parameter_set"]:
        config["parameters"].pop(key)

parameters = config["parameters"][config["parameter_set"]] # short alias for used set of parameters

for tool in parameters["tool_options"]: # sort datatypes in case of mixed datatypes to avoid double calculations
    for option_set in parameters["tool_options"][tool]:
        if "main_datatypes" in parameters["tool_options"][tool][option_set]:
            parameters["tool_options"][tool][option_set]["main_datatypes"] = sorted(parameters["tool_options"][tool][option_set]["main_datatypes"])

for tool in config["other_tool_option_sets"]: # select active set of option for tools other than coretools
    parameters["tool_options"][tool] = parameters["tool_options"][tool][config["other_tool_option_sets"][tool]]

# TODO: check if it is possible to optimize code below

for datatype in config["final_kmer_datatypes"]:
    if (datatype not in config["data_feature_dict"]["fastq"]) and (datatype not in config["data_feature_dict"]["fasta"]):
        raise ValueError(f"ERROR!!! final kmer datatype ({datatype}) is absent among input fastq-based({','.join(config['data_feature_dict']['fastq'])}) and fasta-based({','.join(config['data_feature_dict']['fasta'])}) datatypes")

#check if final_kmer_tool is present in "kmer_counter_list"
if config["final_kmer_counter"] not in parameters["tool_options"]["kmer_qc"]["kmer_counter_list"]:
    parameters["tool_options"]["kmer_qc"]["kmer_counter_list"].append(config["final_kmer_counter"])

#check if final_kmer_length is present in parameters of final_kmer_tool

for dat_type in config["data_feature_dict"]['genome_size']:
    if config["final_kmer_length"] not in parameters["tool_options"][config["final_kmer_counter"]][dat_type]["kmer_length"]:
        parameters["tool_options"][config["final_kmer_counter"]][dat_type]["kmer_length"].append(config["final_kmer_length"])
        #logging.info("Warning! Final_kmer_length is not in parameters of final_kmer_counter! Added...")

for qc_step in "coverage", "merqury", "purge_dups":
    parameters["tool_options"]["assembly_qc"][qc_step]["datatype_list"] = list(set(parameters["tool_options"]["assembly_qc"][qc_step]["datatype_list"]) & set(config["data"].keys()))

#----
#---- Check selected stages ----
brute_logger.info("Preparing pipeline...")
brute_logger.info(TAB + "Selected stages:")
for stage in config["stage_list"]:
    brute_logger.info(TAB * 2  + stage)

for stage in config["stage_list"]:
    if stage not in config["allowed_stage_list"]:
        raise ValueError(f"ERROR!!! Stage {stage} is absent among config['allowed_stage_list']! Check for mistypes... ")

if ("draft_qc" in config["stage_list"]) and ("contig" in config["stage_list"]):
    raise ValueError("ERROR!!! 'draft_qc' and 'contig' are mutually exclusive stage. Select one of them. "
                     "'draft_qc' is a valid choice only if you start from previously created draft assembly. "
                     "If you start from reads, select 'contig'")

stage_dict = OrderedDict()
results_list = []
#---- Configure stages ----

assembler_option_set_group_dict = {} # TODO: this dict is necessary for hifiasm to avoid recalculation of error-corrected reads. Simplify the story. This dict is also used in class Stage and in Hifiasm.smk file

for stage_index in range(0, len(config["stage_list"])):
    stage_dict[config["stage_list"][stage_index]] = Stage(config=config, stage_name=config["stage_list"][stage_index],
                                                          logger=brute_logger,
                                                          prev_stage=None if stage_index == 0 else config["stage_list"][stage_index - 1])
#----
#---- Request files ----
for stage in stage_dict:
    results_list += stage_dict[stage].request_files()
#----

#---- Save configuration and input files ----
final_config_yaml = config["out_dir"] / "config/config.final.yaml"
final_data_yaml = config["out_dir"] / "config/final_data.yaml"

with open(final_config_yaml, 'w') as final_config_fd, open(final_data_yaml, 'w') as final_data_fd:
    yaml.dump(convert_posixpath2str_in_dict(config), final_config_fd, default_flow_style=False, sort_keys=False)
    yaml.dump(convert_posixpath2str_in_dict(config["data"]), final_data_fd, default_flow_style=False, sort_keys=False)

#-------------------------------------------
localrules: all


#---- Global wildcard constrains ----
wildcard_constraints:
    se_datatype="|".join(["hifi", "lqccs", "nanopore", "simplex", "duplex", "ultralongnano", "adaptivenano"]),
    pe_datatype="|".join(["illumina", "hic"]),
    longread_datatype="|".join(["hifi", "lqccs", "nanopore", "simplex", "duplex", "ultralongnano", "adaptivenano"]),
    pacbio_datatype="|".join(["hifi", "lqccs"]),
    nanopore_datatype="|".join(["nanopore", "simplex", "duplex", "ultralongnano", "adaptivenano"]),
    fastqc_datatype="|".join(config["data_feature_dict"]["fastqc"]),
    fileprefix="[^/]+",
    pairprefix="[^/]+",
    fasta_prefix="[^/]+",
    len_prefix="[^/]+",
    scaffold_length="[^./]+",
    datatype="[^./]+", #can be a complex datatype, for example 'hifi_simplex'
    dype="[^_./]+", # elementary datatype only, for example 'hifi' or 'simplex', no '_' is allowed
    haplotype="[^_./@]+",
    merged_haplotype="combined|reordered",
    stage="[^/]+",
    assembly_stage="[^/]+",
    kmer_length="[0-9]+",
    kmer_tool="[^.]+",
    meryl_db=".*meryl.*",
    phasing_kmer_length="[^./]+", # can be an int number or 'NA' in case of no phasing
    genome_prefix="[^/]+",
    correction_options="[^/]+",
    gfa_prefix="[^/]*hap[^/]*|[^/]*alt[^/]*",
    gfa_dir=".*contig.*",
    parameters="[^/]+",
    parameters_prefix="[^/]+",
    prev_stage_parameters="[^/]+",
    purge_dups_parameters="[^/]+",
    dedup_parameters="[^/]+",
    hic_scaffolding_parameters="[^/]+",
    gap_closing_parameters="[^/]+",
    busco5_lineage="[^/]+odb10|[^/]+odb12",
    busco6_lineage="[^/]+odb12\.[^/]*",
    busco_lineage="[^/]*",
    window="[0-9]+",
    step="[0-9]+",
    reference="[^/]+",
    query_prefix="[^/]+",
    target_prefix="[^/]+",
    tab_file_prefix="[^/]+",
    min_target_len="[^0][0-9]+",
    query_scaffold_length="[^/]+",
    target_scaffold_length="[^/]+",
    cov_settings="|".join(parameters["tool_options"]["mosdepth"]["options"].keys()),
    cov_type="[^./]+",
    mapq="[0-9]+",
    min_mapq="[0-9]+",
    resolution="[0-9]+",
    pretext_res="default|low_res|high_res|ultra_res",
    track_type="[^./]+",
    threshold_type="[^/]+"

#---- Final rule ----
pd.Series(results_list).to_csv(config["out_dir"] / "requested_files.tab", sep="\t", header=False, index=False)
rule all:
    input:
        results_list
#----

#---- Include section ----
include: "workflow/rules/General/Log.smk" # DONE
include: "workflow/rules/General/Links.smk" # DONE
include: "workflow/rules/Preprocessing/Files.smk" # DONE
include: "workflow/rules/Preprocessing/Combine.smk" # DONE # TODO: probably not well tested
include: "workflow/rules/Tools/QCFiltering/FastQC.smk" # DONE
include: "workflow/rules/Tools/QCFiltering/MultiQC.smk" # DONE
include: "workflow/rules/Tools/QCFiltering/NanoQC.smk" # DONE
include: "workflow/rules/Tools/QCFiltering/NanoPlot.smk" # DONE
include: "workflow/rules/Tools/QCFiltering/TADbit.smk" # DONE

include: "workflow/rules/Tools/QCFiltering/Cutadapt.smk" # DONE
include: "workflow/rules/Tools/QCFiltering/HiCTrim.smk" # DONE
include: "workflow/rules/Tools/QCFiltering/Trimmomatic.smk" # DONE
include: "workflow/rules/Tools/QCFiltering/Final.smk" # DONE     # TODO: probably not well tested
include: "workflow/rules/Tools/QCFiltering/Nanopore.smk"    # TODO: refactored, but not tested
include: "workflow/rules/Tools/Contamination/Kraken2.smk"  # TODO: refactored, but not tested

include: "workflow/rules/Tools/Kmer/Jellyfish.smk" # DONE
include: "workflow/rules/Tools/Kmer/Meryl.smk"    # DONE
include: "workflow/rules/Tools/Kmer/Yak.smk"      # DONE
include: "workflow/rules/Tools/Kmer/Smudgeplot.smk" # TODO: refactor
include: "workflow/rules/Tools/Kmer/GCplot.smk"     # TODO: refactor
include: "workflow/rules/Tools/Kmer/Genomescope.smk"  # DONE
include: "workflow/rules/Tools/Kmer/Krater.smk"       # DONE

include: "workflow/rules/Stages/contig/Common.smk" # DONE
include: "workflow/rules/Stages/contig/Hifiasm.smk" # DONE
include: "workflow/rules/Stages/contig/Verkko.smk" # TEST
#include: "workflow/rules/Stages/contig/NextDenovo.smk" # DONE
include: "workflow/rules/Stages/contig/Flye.smk" # DONE
include: "workflow/rules/Tools/Graph/GFA.smk" # DONE
include: "workflow/rules/Tools/Contamination/FCS.smk" # TODO: after refactoring, only a case with skip_fcs=True was tested so far

include: "workflow/rules/Tools/General/Sequence.smk" # DONE
include: "workflow/rules/Stats/General.smk" # DONE
include: "workflow/rules/QCAssembly/BUSCO.smk" # # DONE
include: "workflow/rules/Tools/BUSCO/BUSCO5.smk" # DONE
include: "workflow/rules/Tools/BUSCO/BUSCO6.smk" # # DONE
include: "workflow/rules/QCAssembly/QUAST.smk" # DONE
include: "workflow/rules/QCAssembly/Merqury.smk" # DONE
include: "workflow/rules/QCAssembly/GCTrack.smk" # DONE
include: "workflow/rules/QCAssembly/SangerTelomereTrack.smk" # DONE
include: "workflow/rules/Tools/Telomere/SangerTelomere.smk" # DONE
include: "workflow/rules/QCAssembly/TelomereTidkTrack.smk"  # DONE
include: "workflow/rules/Tools/Telomere/Tidk.smk" # DONE   # TODO: not well tested

include: "workflow/rules/QCAssembly/Curation.smk" # DONE     # TODO: probably not well tested
include: "workflow/rules/QCAssembly/GapTrack.smk" # DONE
include: "workflow/rules/QCAssembly/WindowmaskerTrack.smk" # DONE
include: "workflow/rules/Tools/Repeats/Windowmasker.smk" # DONE
include: "workflow/rules/QCAssembly/TRFTrack.smk" # DONE
include: "workflow/rules/Tools/Repeats/TRF.smk" # DONE
include: "workflow/rules/QCAssembly/CoverageTrack.smk" # DONE     # TODO: needs testing on illumina data
include: "workflow/rules/QCAssembly/Purge_dups.smk" # DONE     # TODO: probably not well tested
include: "workflow/rules/Tools/Deduplication/Purge_dups.smk"  # DONE     # TODO: probably not well tested
include: "workflow/rules/Tools/RefScaffolding/RagTag.smk" # DONE

include: "workflow/rules/QCAssembly/CombineHaplotypes.smk" # DONE # TODO: be ready it might interfere with other rules, for example, gfa2fasta
include: "workflow/rules/QCAssembly/MicroChromosomes.smk" # DONE
include: "workflow/rules/QCAssembly/HiCmap.smk" # DONE
include: "workflow/rules/QCAssembly/HiGlass.smk" # DONE
include: "workflow/rules/QCAssembly/Pretext.smk" # DONE # TODO: probably not well tested
include: "workflow/rules/Tools/HiC/Pretext.smk" # DONE # TODO: probably not well tested
include: "workflow/rules/QCAssembly/PretextPerChr.smk" # DONE # TODO: probably not well tested

include: "workflow/rules/Tools/Repeats/Masking.smk" # DONE
include: "workflow/rules/Tools/WGA/LAST.smk"     # DONE

include: "workflow/rules/Stages/read_phasing/ReadPhasing.smk" # DONE

include: "workflow/rules/Tools/Alignment/Index.smk" # DONE
include: "workflow/rules/Tools/Alignment/Common.smk" # DONE # TODO: probably not well tested
include: "workflow/rules/Tools/Alignment/Stats.smk" # DONE
include: "workflow/rules/Tools/Alignment/PostAlignment.smk" # TODO: finish refactoring

if "hic_alignment" in stage_dict:
    include: "workflow/rules/Stages/hic_alignment/Common.smk" # DONE
    include: "workflow/rules/Stages/hic_alignment/Arima.smk" # DONE
    include: "workflow/rules/Stages/hic_alignment/BWA.smk" # DONE
    include: "workflow/rules/Stages/hic_alignment/Pairtools.smk" # DONE
    include: "workflow/rules/Stages/hic_alignment/Juicer.smk" # DONE

if "hic_scaffolding" in stage_dict:
    include: "workflow/rules/Stages/hic_scaffolding/YAHS.smk"
    include: "workflow/rules/Stages/hic_scaffolding/3DDNA.smk" # DONE

include: "workflow/rules/Tools/Deduplication/Hapsolo.smk"     # TODO: test
if "dedup" in stage_dict:
    include: "workflow/rules/Stages/dedup/Common.smk"
    include: "workflow/rules/Stages/dedup/HapSolo.smk"       # TODO: test
    include: "workflow/rules/Stages/dedup/Purge_dups.smk"
    include: "workflow/rules/Stages/dedup/ComboPurge.smk"    # TODO: test

if "ref_scaffolding" in stage_dict:
    pass
    #include: "workflow/rules/ref_scaffolding/RagTag.smk" # TODO: refactor

if "gap_closing" in config["stage_list"]:
    include: "workflow/rules/Stages/gap_closing/Samba.smk" # TODO: test refactored code on big genomes

if "polishing" in stage_dict:
    include: "workflow/rules/Stages/polishing/NextPolish2.smk" # TODO: test

include: "workflow/rules/Tools/Conversion/Bam2bed.smk" # TODO: not tested
include: "workflow/rules/Tools/Alignment/Winnowmap.smk" # TODO: not tested
include: "workflow/rules/Stages/mtdna/MitoHiFi.smk" # DONE
include: "workflow/rules/Stages/mtdna/Mitoz.smk"    # TODO: check code

"""

    include: "workflow/rules/Purge_dups/Purge_dupsQC.smk"
"""