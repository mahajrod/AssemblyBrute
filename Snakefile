import os
import sys

import yaml
#import logging
import shutil
from copy import deepcopy
from collections import OrderedDict
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path, PosixPath

import pandas as pd

#---- Read config files ----
#-------- Read core config file --------
with open(config["main_config_file"], "r") as core_yaml_fd:
    config.update(yaml.safe_load(core_yaml_fd))
#---------------------------------------
#-------- Read resources config files --------
for resource, res_datatype in zip(["threads", "memory_mb", "time"], [int, int, str]):
    resource_df = pd.read_csv(config["resources"][resource], sep="\t", header=0, index_col=0)
    for config_label in resource_df.columns:
        config["parameters"][config_label][resource] = resource_df[config_label].to_dict(OrderedDict)

#---------------------------------------------

#---------------------------

#---- Include sections for functions ----
include: "workflow/functions/option_parsing.py"
include: "workflow/functions/general_parsing.py"
#----------------------------------------


#-- Initialization of path variables from config file --
#logging.info("Initialization of path variables...")
#---- Initialization of path variables for input----
input_dir_path = Path(config["input_dir"])

input_dict = {}
data_types = config["data_types"].split(",")

for datatype in data_types:
    input_dict[datatype] = {}
    input_dict[datatype]["dir"] = input_dir_path / datatype
    input_dict[datatype]["run_dir"] = input_dict[datatype]["dir"] / "run"
    if datatype == "bionano":
        input_dict[datatype]["cmap"] = None # TODO: implement parsing bionano .cmap filename
    else:
        input_dict[datatype]["fastq_dir"] = input_dict[datatype]["dir"] / "fastq"
        input_dict[datatype]["fasta_dir"] = input_dict[datatype]["dir"] / "fasta"
#----
#---- Initialization of path variables for output ----
out_dir_path = Path(config["out_dir"])
output_dict = {}

for first_level_sub_dir in config["first_level_subdir_list"]:
    output_dict[first_level_sub_dir] = out_dir_path / first_level_sub_dir

#----
#---- Initialization path variables for resources ----
#----
#---- Setting mode of pipeline ----
#.info("Setting and adjusting pipeline mode...")

#pipeline_mode = config["mode"]
#starting_point = config["starting_point"]

#-------- Verification of input datatypes --------

fastq_based_data_type_set = set(data_types) & set(config["fastq_based_data"])
fasta_based_data_type_set = set(data_types) & set(config["fasta_based_data"])
fastqc_data_type_set = fastq_based_data_type_set & set(config["fastqc_data_types"])
long_read_data_type_set = set(data_types) & set(config["long_read_data"])
genome_size_estimation_data_type_set = set(config["genome_size_estimation_data"]) & fastq_based_data_type_set & set(data_types)
coverage_track_data_type_set = set(data_types) & set(config["coverage_track_data"])
variant_calling_data_type_set = set(data_types) & set(config["variant_calling_data"])

#logging.info("Verifying datatypes...")
for d_type in data_types:
    if d_type not in config["allowed_data_types"]:
        #logging.error("Unknown data type: {0}".format(d_type))
        raise ValueError("ERROR!!! Unknown data type: {0}".format(d_type))

if config["final_kmer_datatype"] not in fastq_based_data_type_set:
    raise ValueError("ERROR!!! final_kmer_datatype ({0}) is absent among input fastq-based datatypes({1})".format(config["final_kmer_datatype"],
                                                                                                                  ",".join(fastq_based_data_type_set)))

#--------

#----

#---- Checking input files ----
#logging.info("Checking input files...")

input_filedict = {}
input_file_prefix_dict = {}
input_fasta_filedict = {}
input_fasta_file_prefix_dict = {}
input_reference_filedict = {}

input_forward_suffix_dict = {}
input_reverse_suffix_dict = {}
input_pairprefix_dict = {}

for d_type in fastq_based_data_type_set:
    input_filedict[d_type] = find_fastqs(input_dict[d_type]["fastq_dir"], fastq_extension=config["fastq_extension"])
    input_file_prefix_dict[d_type] = list(map(lambda s: str(s.name)[:-len(config["fastq_extension"])],
                                                input_filedict[d_type]))

for d_type in fasta_based_data_type_set:
    input_fasta_filedict[d_type] = find_fastas(input_dict[d_type]["fasta_dir"], fasta_extension=config["fasta_extension"])
    input_fasta_file_prefix_dict[d_type] = list(map(lambda s: str(s.name)[:-len(config["fasta_extension"])],
                                                input_fasta_filedict[d_type]))

#---- detect datatypes and check if datatype has files in both fasta and fastq formats ----
datatype_format_dict = {}
datatype_extension_dict = {}
for d_type in set(data_types):
    if (d_type in fastq_based_data_type_set) and (d_type in fasta_based_data_type_set):
        if (len(input_fasta_filedict[d_type]) > 0) and (len(input_filedict[d_type]) > 0):
            raise  ValueError("Error!!! Datatype {0} has input files in both fastq ({1}) and fasta ({2}) formats!".format(d_type,
                                                                                                                          " ".join(input_filedict[d_type]),
                                                                                                                          " ".join(input_fasta_filedict[d_type])))
        elif len(input_fasta_filedict[d_type]) > 0:
            datatype_format_dict[d_type] = "fasta"
            datatype_extension_dict[d_type] = config["fasta_extension"]
        elif len(input_filedict[d_type]) > 0:
             datatype_format_dict[d_type] = "fastq"
             datatype_extension_dict[d_type] = config["fastq_extension"]

if "reference" in set(data_types):
    reference_input_dir = input_dict[datatype]["dir"]
    reference_genomes_list = [element.name for element in reference_input_dir.glob("*")]
    for genome in reference_genomes_list:
        input_reference_filedict[genome] = {}
        for filetype in "fasta", "syn", "whitelist", "orderlist":
            input_reference_filedict[genome][filetype] = list((reference_input_dir / genome).glob("*.{0}".format(filetype)))
            if len(input_reference_filedict[genome][filetype]) > 1:
                raise ValueError("ERROR!!! There is more than one {0} file for reference {1}".format(filetype, genome))
            input_reference_filedict[genome][filetype] = input_reference_filedict[genome][filetype][0]
        #for filetype in "mtdna.fasta", "mtdna.gb":
        #    input_reference_filedict[genome][filetype] = list((reference_input_dir / genome / "mtdna").glob("*.{0}".format(filetype)))
        #    if len(input_reference_filedict[genome][filetype]) > 1:
        #        raise ValueError("ERROR!!! There is more than one {0} file for reference {1}".format(filetype, genome))
        #    input_reference_filedict[genome][filetype] = input_reference_filedict[genome][filetype][0]
    #print(input_reference_filedict)
#------------------------------------------------------------------------------------------

# check filenames of paired data
for d_type in set(config["paired_fastq_based_data"]) & fastq_based_data_type_set:
   if (len(input_filedict[d_type]) % 2) != 0:
        raise ValueError("ERROR!!! {0} fastq files seems to be unpaired or misrecognized".format(d_type))
   for forward, reverse in zip(input_filedict[d_type][::2], input_filedict[d_type][1::2]):
        if p_distance(str(forward), str(reverse), len(str(forward))) > 1:
            raise ValueError("ERROR!!! Forward and reverse read files differs by more than one symbol:\n\t{0}\n\t{1}".format(str(forward),
                                                                                                                             str(reverse)))
#get_suffixes for paired fastq data
for d_type in set(config["paired_fastq_based_data"]) & fastq_based_data_type_set:
    input_forward_suffix_dict[d_type] = set()
    input_reverse_suffix_dict[d_type] = set()
    input_pairprefix_dict[d_type] = []
    for forward_prefix, reverse_prefix in zip(input_file_prefix_dict[d_type][::2], input_file_prefix_dict[d_type][1::2]):
        common_prefix, forward_suffix, reverse_suffix = get_common_prefix_ans_suffixes(forward_prefix, reverse_prefix)
        input_pairprefix_dict[d_type].append(common_prefix)
        input_forward_suffix_dict[d_type].add(forward_suffix)
        input_reverse_suffix_dict[d_type].add(reverse_suffix)
    if (len(input_forward_suffix_dict[d_type]) > 1) or (len(input_reverse_suffix_dict[d_type]) > 1):
        raise ValueError("ERROR!!! Multiple different suffixes in filenames of %s data!" % d_type)

    input_forward_suffix_dict[d_type] = list(input_forward_suffix_dict[d_type])[0]
    input_reverse_suffix_dict[d_type] = list(input_reverse_suffix_dict[d_type])[0]

#for d_type in fastq_based_data_type_set: # add prefixes of files for se data to simplify dealing with wildcards
#    if d_type not in input_pairprefix_dict:
#        input_pairprefix_dict[d_type] = input_file_prefix_dict[d_type]

if "bionano" in data_types: # TODO: modify when input for bionano will be clear
    input_filedict["bionano"] = find_cmap(input_dict["bionano"]["dir"], cmap_extension=config["cmap_extension"])


#---- Initialize tool parameters ----
#logging.info("Initializing tool parameters...")
#check if custom restriction sites were provided:
if config["custom_enzyme_set"] is not None:
    config["parameters"]["default"]["tool_options"]["salsa2"]["restriction_seq"]["custom"] = config["custom_enzyme_set"]
    if "tool_options" in config["parameters"][config["parameter_set"]]:
        if "salsa2" in config["parameters"][config["parameter_set"]]["tool_options"]:
            if "restriction_seq" in config["parameters"][config["parameter_set"]]["tool_options"]["salsa2"]:
                if "custom" in config["parameters"][config["parameter_set"]]["tool_options"]["salsa2"]["restriction_seq"]:
                    config["parameters"][config["parameter_set"]]["tool_options"]["salsa2"]["restriction_seq"]["custom"] = config["custom_enzyme_set"]
    config["hic_enzyme_set"] = "custom"

    if config["custom_enzyme_set_is_no_motif"]: # register custom enzyme as producicing no ligation motives
        config["no_motif_enzyme_sets"].append(custom)

if config["parameter_set"] not in config["parameters"]:
    raise ValueError("Error!!! Unknown set of tool parameters: {0}".format(config["parameter_set"]))

copy_absent_entries(config["parameters"]["default"], config["parameters"][config["parameter_set"]]) # set default values for options absent in  "parameter_set"

for key in list(config["parameters"].keys()): # remove unused sets of parameters
    if key != config["parameter_set"]:
        config["parameters"].pop(key)

parameters = config["parameters"][config["parameter_set"]] # short alias for used set of parameters

for tool in config["other_tool_option_sets"]: # select active set of option for tools other than coretools
    parameters["tool_options"][tool] = parameters["tool_options"][tool][config["other_tool_option_sets"][tool]]

#check if final_kmer_tool is present in "kmer_counter_list"
if config["final_kmer_counter"] not in config["kmer_counter_list"]:
    config["kmer_counter_list"].append(config["final_kmer_counter"])
    #logging.info("Warning! final_kmer_counter is not in kmer_counter_list! Added...")

#check if final_kmer_length is present in parameters of final_kmer_tool
for dat_type in genome_size_estimation_data_type_set:
    if config["final_kmer_length"] not in parameters["tool_options"][config["final_kmer_counter"]][dat_type]["kmer_length"]:
        parameters["tool_options"][config["final_kmer_counter"]][dat_type]["kmer_length"].append(config["final_kmer_length"])
        #logging.info("Warning! Final_kmer_length is not in parameters of final_kmer_counter! Added...")

#Kraken scan datatype
kraken_scan_data_type_set = set(data_types) & set(config["kraken_scan_data"])

#----
#---- Configure stages ----
config["stage_list"] = []

# Select configuration and combine stages from all mega_stages in a single list without nesting
if config["mode"] == "preprocessing":
    mega_stage_list = ["preprocessing"]
elif config["mode"] == "qc":
    mega_stage_list = ["preprocessing", "qc"]
elif config["mode"] == "assembly":
    mega_stage_list = ["preprocessing", "qc", "assembly"]
elif config["mode"] == "finalization":
    mega_stage_list = ["preprocessing", "qc", "finalization"]
else:
    raise ValueError("ERROR!!! Unknown mode: %s" % config["mode"])

for mega_stage in mega_stage_list:
    custom_megastage_entry = "custom_" + mega_stage + "_stages"
    if (custom_megastage_entry in config) and (config[custom_megastage_entry]):
        config["stage_list"].append(config[custom_megastage_entry])
    else:
        config["stage_list"] += config["allowed_stage_list"][mega_stage][config[mega_stage + "_mode"]][config["starting_point"]]

stage_dict = OrderedDict()
for stage, stage_index in zip(config["stage_list"], range(0, len(config["stage_list"]))):
    stage_dict[stage] = OrderedDict()
    stage_dict[stage]["prev_stage"] = None if stage_index == 0 else config["stage_list"][stage_index-1]

#----
#---- Save configuration and input files ----
final_config_yaml = output_dict["config"] / "config.final.yaml"
final_input_yaml = output_dict["config"] / "input.final.yaml"

os.makedirs(output_dict["config"], exist_ok=True)

with open(final_config_yaml, 'w') as final_config_fd, open(final_input_yaml, 'w') as final_input_fd:
    yaml.dump(convert_posixpath2str_in_dict(config), final_config_fd, default_flow_style=False, sort_keys=False)
    yaml.dump(convert_posixpath2str_in_dict(input_dict), final_input_fd, default_flow_style=False, sort_keys=False)

#-------------------------------------------
localrules: all
#ruleorder: create_fastq_links > fastqc

results_dict = {}

haplotype_list = ["hap{0}".format(i) for i in range(1, config["ploidy"] + 1)] # TODO: obsolete: remove and fix issues
primary_haplotype = "hap1" # TODO: obsolete: remove and fix issues

results_list = []

for conda_env in config["conda"]:
    if "pip" in config["conda"][conda_env]:
        results_list += [expand("results/config/pip.{conda_env}.requirements",
                                conda_env=[conda_env])]

#---- Create output filelist ----
if "check_reads" in config["stage_list"]:
    results_list += [
                     final_config_yaml,
                     final_input_yaml
                     ]

if "check_draft" in config["stage_list"]:
    results_list += [ ] # TODO: implement


if ("read_qc" in config["stage_list"]) and (not config["skip_read_qc"]):
    results_list += [*[expand(output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip",
                               datatype=[dat_type, ],
                               stage=["raw", ],
                               fileprefix=input_file_prefix_dict[dat_type],) for dat_type in fastqc_data_type_set ],
                      expand(output_dict["qc"] / "multiqc/{datatype}/{stage}/multiqc.{datatype}.{stage}.report.html",
                             datatype=fastqc_data_type_set ,
                             stage=["raw",]),
                      *[expand(output_dict["qc"] / "nanoplot/{datatype}/{stage}/{fileprefix}.Yield_By_Length.png",
                               datatype=[dat_type, ],
                               stage=["raw", ],
                               fileprefix=input_file_prefix_dict[dat_type],) for dat_type in long_read_data_type_set],
                    *[expand(output_dict["qc"] / "nanoqc/{datatype}/{stage}/{fileprefix}",
                               datatype=[dat_type, ],
                               stage=["raw", ],
                               fileprefix=input_file_prefix_dict[dat_type],) for dat_type in long_read_data_type_set],
                     ]


if "draft_qc" in config["stage_list"]:
    draft_file_dict = get_input_assemblies(input_dir_path / "draft/fasta", config["ploidy"], config["assembly_fasta_extension"])
    stage_dict["draft_qc"]["parameters"] = {}

    for qcer in config["stage_coretools"]["draft_qc"]["default"]:
        for option_set in config["coretool_option_sets"][qcer]:
            parameters_label="{0}_{1}".format(qcer, option_set)
            stage_dict["draft_qc"]["parameters"][parameters_label] = {}
            stage_dict["draft_qc"]["parameters"][parameters_label]["qcer"] = qcer
            stage_dict["draft_qc"]["parameters"][parameters_label]["option_set"] = {}
            stage_dict["draft_qc"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] = config["ploidy"]
            stage_dict["draft_qc"]["parameters"][parameters_label]["haplotype_list"] = ["hap{0}".format(i) for i in range(1, stage_dict["draft_qc"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] + 1)] if stage_dict["draft_qc"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] > 1 else ["hap0"]
            stage_dict["draft_qc"]["parameters"][parameters_label]["option_set_group"] = None

    parameters_list = list(stage_dict["draft_qc"]["parameters"].keys())

    results_list += [expand(out_dir_path / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["draft_qc"],),
                     *[expand(out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
                                assembly_stage=["draft_qc"],
                                parameters=[parameters_label],
                                genome_prefix=[config["genome_prefix"], ],
                                haplotype=stage_dict["draft_qc"]["parameters"][parameters_label]["haplotype_list"]
                                ) for parameters_label in parameters_list],
                     ]
    if not config["skip_busco"]:
        results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.busco5.{busco_lineage}.tar.gz",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["draft_qc"],
                                haplotype=stage_dict["draft_qc"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["draft_qc"],
                                haplotype=stage_dict["draft_qc"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         ]
    #TODO: remove after debugging
    """
    results_list += [ *[expand(out_dir_path / "{stage}/{parameters}/kmer/{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}",
                               stage=["draft_qc"],
                              parameters=[parameters_label],
                              genome_prefix=[config["genome_prefix"], ],
                              haplotype=stage_dict["draft_qc"]["parameters"][parameters_label]["haplotype_list"],
                              assembly_kmer_length=[31]) for parameters_label in parameters_list],
                      *[expand(out_dir_path / "{stage}/{parameters}/fasta/{haplotype}/{assembly_kmer_length}/{datatype}/{fileprefix}.fasta.gz",
                                stage=["draft_qc"],
                              parameters=[parameters_label],
                                datatype=[config["gap_closing_datatype"]],
                                fileprefix=input_file_prefix_dict[config["gap_closing_datatype"]] if datatype_format_dict[config["gap_closing_datatype"]] == "fastq" else input_fasta_file_prefix_dict[config["gap_closing_datatype"]],
                              genome_prefix=[config["genome_prefix"], ],
                              haplotype=stage_dict["draft_qc"]["parameters"][parameters_label]["haplotype_list"],
                              assembly_kmer_length=[31]) for parameters_label in parameters_list],

                      ]
    """
    if "gap_closing" in config["stage_list"]:
        prev_stage = "draft_qc"

        gap_closer_list = config["stage_coretools"]["gap_closing"]["default"]
        stage_dict["gap_closing"]["parameters"] = {}

        for gap_closer in gap_closer_list:
            for option_set in config["coretool_option_sets"][gap_closer]:
                for prev_parameters in stage_dict[prev_stage]["parameters"]:
                    parameters_label = "{0}..{1}_{2}".format(prev_parameters, gap_closer, option_set)
                    stage_dict["gap_closing"]["parameters"][parameters_label] = {}
                    stage_dict["gap_closing"]["parameters"][parameters_label]["included"] = True
                    stage_dict["gap_closing"]["parameters"][parameters_label]["gap_closer"] = gap_closer
                    stage_dict["gap_closing"]["parameters"][parameters_label]["prev_stage"] = prev_stage
                    stage_dict["gap_closing"]["parameters"][parameters_label]["option_set"] = deepcopy(parameters["tool_options"][gap_closer][option_set])
                    stage_dict["gap_closing"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] = config["ploidy"]
                    stage_dict["gap_closing"]["parameters"][parameters_label]["haplotype_list"] = ["hap{0}".format(i) for i in range(1, stage_dict["gap_closing"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] + 1)] if stage_dict["gap_closing"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] > 1 else ["hap0"]
                    stage_dict["gap_closing"]["parameters"][parameters_label]["option_set_group"] = None

        parameters_list = list(stage_dict["gap_closing"]["parameters"].keys())
        results_list += [*[expand(out_dir_path / "gap_closing/{parameters}/{genome_prefix}.gap_closing.{haplotype}.len",
                                parameters=[parameters_label],
                                genome_prefix=[config["genome_prefix"], ],
                                haplotype=stage_dict["gap_closing"]["parameters"][parameters_label]["haplotype_list"]
                                ) for parameters_label in parameters_list],
                         expand(out_dir_path / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["gap_closing"],),
                         [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.pre.hic",
                                  assembly_stage=["gap_closing"],
                                  parameters=[parameters_label],
                                  genome_prefix=[config["genome_prefix"], ],
                                  haplotype=stage_dict["gap_closing"]["parameters"][parameters_label]["haplotype_list"],
                                  phasing_kmer_length=[stage_dict["gap_closing"]["parameters"][parameters_label]["option_set"]["phasing_kmer_length"]])
                           for parameters_label in parameters_list] if not config["skip_hic_file"] else []
                         ]

        if not config["skip_busco"]:
            results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.busco5.{busco_lineage}.tar.gz",
                                    busco_lineage=config["busco_lineage_list"],
                                    genome_prefix=[config["genome_prefix"], ],
                                    assembly_stage=["gap_closing"],
                                    haplotype=stage_dict["gap_closing"]["parameters"][parameters_label]["haplotype_list"],
                                    parameters=[parameters_label]) for parameters_label in parameters_list],
                             *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.busco.merged.tsv",
                                    busco_lineage=config["busco_lineage_list"],
                                    genome_prefix=[config["genome_prefix"], ],
                                    assembly_stage=["gap_closing"],
                                    haplotype=stage_dict["gap_closing"]["parameters"][parameters_label]["haplotype_list"],
                                    parameters=[parameters_label]) for parameters_label in parameters_list],
                         ]

if ("filter_reads" in config["stage_list"]) and (not config["skip_filter_reads"]):
    results_list += [expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
                            fileprefix=input_file_prefix_dict["hifi"]) if "hifi" in fastq_based_data_type_set else [],
                    expand(output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip",
                           datatype=["hifi", ],
                           stage=["filtered", ],
                           fileprefix=input_file_prefix_dict["hifi"],
                           ) if "hifi" in fastq_based_data_type_set else [],
                    expand(output_dict["qc"] / "multiqc/{datatype}/{stage}/multiqc.{datatype}.{stage}.report.html",
                           datatype=["hifi"],
                           stage=["filtered",]) if "hifi" in fastq_based_data_type_set else [],
                    *[[expand(output_dict["kmer"] / "{datatype}/{stage}/genomescope/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.genomescope.parameters",
                           datatype=[dat_type,],
                           genome_prefix=[config["genome_prefix"], ],
                           stage=["filtered",],
                           kmer_tool=[kmer_tool,],
                           kmer_length=parameters["tool_options"][kmer_tool][dat_type]["kmer_length"],
                           ) for kmer_tool in config["kmer_counter_list"] ]  for dat_type in genome_size_estimation_data_type_set],
                    ]
    if not config["skip_nanoqc"]:
        results_list += [

                        *[expand(output_dict["qc"] / "nanoqc/{datatype}/{stage}/{fileprefix}",
                                   datatype=[dat_type, ],
                                   stage=["filtered", ],
                                   fileprefix=input_file_prefix_dict[dat_type],) for dat_type in long_read_data_type_set],
                        expand(output_dict["qc"] / "nanoqc/{datatype}/{stage}/{fileprefix}",
                                   datatype=["nanopore", ],
                                   stage=["trimmed", ],
                                   fileprefix=input_file_prefix_dict["nanopore"],) if "nanopore" in long_read_data_type_set else [],
                        ]
    if not config["skip_nanoplot"]:
        results_list += [*[expand(output_dict["qc"] / "nanoplot/{datatype}/{stage}/{fileprefix}.Yield_By_Length.png",
                               datatype=[dat_type, ],
                               stage=["filtered", ],
                               fileprefix=input_file_prefix_dict[dat_type],) for dat_type in long_read_data_type_set],
                        expand(output_dict["qc"] / "nanoplot/{datatype}/{stage}/{fileprefix}.Yield_By_Length.png",
                                   datatype=["nanopore", ],
                                   stage=["trimmed", ],
                                   fileprefix=input_file_prefix_dict["nanopore"],) if "nanopore" in long_read_data_type_set else [],
                        ]

    if config["database_set"]["kraken2"] and kraken_scan_data_type_set and (not config["skip_kraken"]):
        results_list += [expand(out_dir_path / "contamination_scan/kraken2/{datatype}/kraken2.{database}.report",
                               datatype=kraken_scan_data_type_set,
                               database=config["database_set"]["kraken2"],
                               )
                        ]

if "smudgeplot" in config["stage_list"]:
    results_list += [*[[expand(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}_warnings.txt",
                               lower_boundary=parameters["tool_options"]["smudgeplot"]["lower_boundary"],
                               upper_boundary=parameters["tool_options"]["smudgeplot"]["upper_boundary"],
                               datatype=[dat_type,],
                               stage=["filtered",],
                               kmer_tool=[kmer_tool,],
                               kmer_length=parameters["tool_options"][kmer_tool][dat_type]["kmer_length"],
                               ) for kmer_tool in config["kmer_counter_list"] ]  for dat_type in genome_size_estimation_data_type_set],
                    #*[[expand(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.L{lower_boundary}.U{upper_boundary}.subset.kmer.gz",
                    #          lower_boundary=parameters["tool_options"]["smudgeplot"]["lower_boundary"],
                    #          upper_boundary=parameters["tool_options"]["smudgeplot"]["upper_boundary"],
                    #          datatype=[dat_type,],
                    #          stage=["filtered",],
                    #          kmer_tool=[kmer_tool,],
                    #          kmer_length=parameters["tool_options"][kmer_tool][dat_type]["kmer_length"],
                    #          ) for kmer_tool in config["kmer_counter_list"] ]  for dat_type in genome_size_estimation_data_type_set],
                    *[[expand(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.smudgeplot.boundaries",
                              datatype=[dat_type,],
                              stage=["filtered",],
                              kmer_tool=[kmer_tool,],
                              kmer_length=parameters["tool_options"][kmer_tool][dat_type]["kmer_length"],
                              ) for kmer_tool in config["kmer_counter_list"] ]  for dat_type in genome_size_estimation_data_type_set]
                    ]
#if "kat" in config["stage_list"]:
#    results_list += [expand(output_dict["kmer"] / "{datatype}/{stage}/kat/{datatype}.{stage}.{kmer_length}.jellyfish.kat.gcp.mx.png",
#                     datatype=[dat_type,],
#                     stage=["filtered",],
#                     kmer_length=parameters["tool_options"]["kat"][dat_type]["kmer_length"],
#                     )  for dat_type in set(parameters["tool_options"]["kat"]) & set(data_types)
#                    ]
if ("gcp" in config["stage_list"]) and (not config["skip_gcp"]):
    results_list += [expand(output_dict["kmer"] / "{datatype}/{stage}/gcp/{datatype}.{stage}.{kmer_length}.L{min_coverage}.heatmap.png",
                     datatype=[dat_type,],
                     stage=["filtered",],
                     kmer_length=parameters["tool_options"]["gcp"][dat_type]["kmer_length"],
                     min_coverage=parameters["tool_options"]["gcp"][dat_type]["min_coverage"],
                     )  for dat_type in set(parameters["tool_options"]["gcp"]) & set(data_types)
                    ]

if "filter_draft" in config["stage_list"]:
    results_list += [ ] # TODO: implement

if "contig" in config["stage_list"] or "draft_qc" in config["stage_list"]:
    stage_dict["contig"] = {}
    assembler_list = config["stage_coretools"]["contig"][config["contig_datatype"]]
    stage_dict["contig"]["parameters"] = {}
    assembler_option_set_group_dict = {}

    for assembler in assembler_list:
        option_set_group_dict, option_set_group_assignment_dict = None, None
        if assembler == "hifiasm":
            option_set_group_dict, option_set_group_assignment_dict = group_option_sets(parameters["tool_options"]["hifiasm"],
                                                                                        config["tool_specific_features"]["hifiasm"]['options_affecting_error_correction'])
            assembler_option_set_group_dict[assembler] = option_set_group_dict
        for option_set in config["coretool_option_sets"][assembler]:
            parameters_label="{0}_{1}".format(assembler, option_set)
            stage_dict["contig"]["parameters"][parameters_label] = {}
            stage_dict["contig"]["parameters"][parameters_label]["included"] = True
            stage_dict["contig"]["parameters"][parameters_label]["assembler"] = assembler
            stage_dict["contig"]["parameters"][parameters_label]["option_set"] = deepcopy(parameters["tool_options"][assembler][option_set])
            if stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] is None:
               stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] = config["ploidy"]

            stage_dict["contig"]["parameters"][parameters_label]["haplotype_list"] = ["hap{0}".format(i) for i in range(1, stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] + 1)] if stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] > 1 else ["hap0"]
            stage_dict["contig"]["parameters"][parameters_label]["option_set_group"] = option_set_group_assignment_dict[option_set] if option_set_group_assignment_dict is not None else none

            #for option_supergroup in ["options_affecting_error_correction"]:
            #    stage_dict["contig"]["parameters"][parameters_label][option_supergroup] = option_cluster_reverse_dict[assembler][option_supergroup][option_set]
#print (stage_dict)
if "contig" in config["stage_list"]:
    parameters_list = list(stage_dict["contig"]["parameters"].keys())
    #if "hifiasm" in assembler_list:
    #    results_list += [expand(output_dict["error_correction"] / "hifiasm_{correction_options}/{genome_prefix}.contig.ec.bin",
    #                            genome_prefix=[config["genome_prefix"],],
    #                            correction_options=assembler_option_set_group_dict[assembler]),
    #
    #                     ]
    results_list += [*[expand(output_dict["contig"] / "{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
                            genome_prefix=[config["genome_prefix"],],
                            assembly_stage=["contig",],
                            haplotype=stage_dict["contig"]["parameters"][parameters_label]["haplotype_list"] + (["alt" if stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] > 1 else "alt0"] if stage_dict["contig"]["parameters"][parameters_label]["assembler"] == "hifiasm" else []), # TODO: modify "alt" when assemblers other than hifiasm will be added
                            parameters=[parameters_label]) for parameters_label in parameters_list],
                     *[expand(output_dict["contig"] / "{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.unfiltered.gfa.cov",
                            genome_prefix=[config["genome_prefix"],],
                            assembly_stage=["contig",],
                            haplotype=stage_dict["contig"]["parameters"][parameters_label]["haplotype_list"] +  (["alt" if stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] > 1 else "alt0"] if stage_dict["contig"]["parameters"][parameters_label]["assembler"] == "hifiasm" else []), # TODO: modify "alt" when assemblers other than hifiasm will be added
                            parameters=[parameters_label])  for parameters_label in parameters_list],
                     *[expand(output_dict["contig"] / "{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.unfiltered.gfa.lencov",
                            genome_prefix=[config["genome_prefix"],],
                            assembly_stage=["contig",],
                            haplotype=stage_dict["contig"]["parameters"][parameters_label]["haplotype_list"] + (["alt" if stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] > 1 else "alt0"] if stage_dict["contig"]["parameters"][parameters_label]["assembler"] == "hifiasm" else []), # TODO: modify "alt" when assemblers other than hifiasm will be added
                            parameters=[parameters_label]) for parameters_label in parameters_list],
                     *[expand(out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["contig"],
                           haplotype=stage_dict["contig"]["parameters"][parameters_label]["haplotype_list"],
                           parameters=[parameters_label]) for parameters_label in parameters_list],
                    expand(out_dir_path / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["contig"],),
                    ] # Tested only on hifiasm

    if not config["skip_busco"]:
        results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.busco5.{busco_lineage}.tar.gz",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["contig"],
                                haplotype=stage_dict["contig"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["contig"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["contig"],
                                haplotype=stage_dict["contig"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["contig"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         ]
    if (config["tax_id"] is None) or (not config["tax_id"]):
        sys.stderr.write("Tax id was not set, skipping contamination scan in FCS databases...\n")
    else:
        if config["database_set"]["fcs_adaptor"] and (not config["skip_fcs_adaptor"]):
            results_list += [
                            *[expand(out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs_adaptor/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.lenfiltered.{database}.report",
                                   genome_prefix=[config["genome_prefix"], ],
                                   assembly_stage=["contig"],
                                   haplotype=stage_dict["contig"]["parameters"][parameters_label]["haplotype_list"],
                                   parameters=[parameters_label],
                                   database=config["database_set"]["fcs_adaptor"]) for parameters_label in parameters_list],
                            ]
        if config["database_set"]["fcs"] and (not config["skip_fcs"]):
            results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.lenfiltered.{database}.taxonomy",
                                    genome_prefix=[config["genome_prefix"], ],
                                    assembly_stage=["contig"],
                                    haplotype=stage_dict["contig"]["parameters"][parameters_label]["haplotype_list"] + (["alt" if stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] > 1 else "alt0"] if stage_dict["contig"]["parameters"][parameters_label]["assembler"] == "hifiasm" else []), # TODO: modify "alt" when assemblers other than hifiasm will be added
                                    parameters=[parameters_label],
                                    database=config["database_set"]["fcs"]) for parameters_label in parameters_list]
                            ]

if "purge_dups" in config["stage_list"]:
    prev_stage = stage_dict["purge_dups"]["prev_stage"]
    purge_dupser_list = config["stage_coretools"]["purge_dups"]["default"]
    stage_dict["purge_dups"]["parameters"] = {}

    for purge_dupser in purge_dupser_list:
        for option_set in config["coretool_option_sets"][purge_dupser]:
            for prev_parameters in stage_dict[prev_stage]["parameters"]:
                parameters_label = "{0}..{1}_{2}".format(prev_parameters, purge_dupser, option_set)
                stage_dict["purge_dups"]["parameters"][parameters_label] = {}
                stage_dict["purge_dups"]["parameters"][parameters_label]["included"] = True
                stage_dict["purge_dups"]["parameters"][parameters_label]["prev_stage"] = prev_stage
                stage_dict["purge_dups"]["parameters"][parameters_label]["prev_parameters"] = prev_parameters
                stage_dict["purge_dups"]["parameters"][parameters_label]["purge_dupser"] = purge_dupser
                stage_dict["purge_dups"]["parameters"][parameters_label]["option_set"] = parameters["tool_options"][purge_dupser][option_set]
                stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"] = stage_dict[stage_dict["purge_dups"]["prev_stage"]]["parameters"][prev_parameters]["haplotype_list"]

    parameters_list = list(stage_dict["purge_dups"]["parameters"].keys())
    results_list += [
                     *[expand(out_dir_path / "purge_dups/{parameters}/{genome_prefix}.purge_dups.{haplotype}.fasta",
                              genome_prefix=[config["genome_prefix"], ],
                              assembly_stage=["contig"],
                              haplotype=stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"],
                              parameters=[parameters_label]) for parameters_label in parameters_list],
                    *[expand(out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
                             genome_prefix=[config["genome_prefix"], ],
                             assembly_stage=["purge_dups"],
                             haplotype=stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"],
                             parameters=[parameters_label]) for parameters_label in parameters_list],
                    *[expand(out_dir_path /  "{assembly_stage}/{parameters}/assembly_qc/purge_dups/{haplotype}/PB.stat",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["purge_dups"],
                           haplotype=stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"],
                           parameters=[parameters_label]) for parameters_label in parameters_list],
                    expand(out_dir_path /  "{assembly_stage}/{parameters}/assembly_qc/purge_dups/before.comparison.coverage.png",
                           assembly_stage=["purge_dups"],
                           parameters=parameters_list
                           ),
                    expand(out_dir_path / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["purge_dups"],),
                    [[expand(out_dir_path  / "purge_dups/{parameters}/{purge_stage}/{haplotype}/{genome_prefix}.dups.{artefact}.fasta",
                           purge_stage=["first_stage",] if haplotype == "hap0" else ["first_stage", "second_stage"],
                           genome_prefix=[config["genome_prefix"], ],
                           artefact=["junk", "repeat", "haplotig", "ovlp", "highcov"],
                           haplotype=[haplotype],
                           parameters=[parameters_label]) for haplotype in stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"]] for parameters_label in parameters_list],
                    ]
    if not config["skip_busco"]:
        results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.busco5.{busco_lineage}.tar.gz",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["purge_dups", ],
                                haplotype=stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["purge_dups"],
                                #haplotype=stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["purge_dups"],
                                haplotype=stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["purge_dups"],
                                parameters=parameters_list
                                ),
                         ]

if (config["phasing_stage"] in config["stage_list"]) and (not config["skip_phasing"]):

    for datatype in set(data_types) & set(config["read_phasing_data"]):
        if datatype in config["paired_fastq_based_data"]:
            results_list += [*[(expand(out_dir_path / "{stage}/{parameters}/fastq/{haplotype}/{assembly_kmer_length}/{datatype}/{pairprefix}_1.fastq.gz",
                                    datatype=[datatype],
                                    stage=[config["phasing_stage"], ],
                                    parameters=[parameters_label],
                                    pairprefix=input_pairprefix_dict[datatype],
                                    genome_prefix=[config["genome_prefix"], ],
                                    haplotype=stage_dict[config["phasing_stage"]]["parameters"][parameters_label]["haplotype_list"],
                                    assembly_kmer_length=config["assembly_kmer_length"]
                                    ) if len(stage_dict[config["phasing_stage"]]["parameters"][parameters_label]["haplotype_list"]) > 1 else []) for parameters_label in list(stage_dict[config["phasing_stage"]]["parameters"].keys())] ,
                            ]
        else:
            results_list += [*[(expand(out_dir_path / "{stage}/{parameters}/fastq/{haplotype}/{assembly_kmer_length}/{datatype}/{fileprefix}.fastq.gz",
                                    datatype=[datatype],
                                    stage=[config["phasing_stage"], ],
                                    parameters=[parameters_label],
                                    fileprefix=input_file_prefix_dict[datatype],
                                    genome_prefix=[config["genome_prefix"], ],
                                    haplotype=stage_dict[config["phasing_stage"]]["parameters"][parameters_label]["haplotype_list"],
                                    assembly_kmer_length=config["assembly_kmer_length"]
                                    ) if len(stage_dict[config["phasing_stage"]]["parameters"][parameters_label]["haplotype_list"]) > 1 else []) for parameters_label in list(stage_dict[config["phasing_stage"]]["parameters"].keys())],
                            ]
if "hic_scaffolding" in config["stage_list"]:
    prev_stage = stage_dict["hic_scaffolding"]["prev_stage"]
    hic_scaffolder_list = config["stage_coretools"]["hic_scaffolding"]["default"]
    stage_dict["hic_scaffolding"]["parameters"] = {}

    for hic_scaffolder in hic_scaffolder_list:
        for option_set in config["coretool_option_sets"][hic_scaffolder]:
            for prev_parameters in stage_dict[prev_stage]["parameters"]:
                parameters_label = "{0}..{1}_{2}".format(prev_parameters, hic_scaffolder, option_set)
                stage_dict["hic_scaffolding"]["parameters"][parameters_label] = {}
                stage_dict["hic_scaffolding"]["parameters"][parameters_label]["included"] = True
                stage_dict["hic_scaffolding"]["parameters"][parameters_label]["prev_stage"] = prev_stage
                stage_dict["hic_scaffolding"]["parameters"][parameters_label]["prev_parameters"] = prev_parameters
                stage_dict["hic_scaffolding"]["parameters"][parameters_label]["hic_scaffolder"] = hic_scaffolder
                stage_dict["hic_scaffolding"]["parameters"][parameters_label]["option_set"] = parameters["tool_options"][hic_scaffolder][option_set]
                stage_dict["hic_scaffolding"]["parameters"][parameters_label]["haplotype_list"] = stage_dict[stage_dict["hic_scaffolding"]["prev_stage"]]["parameters"][prev_parameters]["haplotype_list"]

                if (len(stage_dict["hic_scaffolding"]["parameters"][parameters_label]["haplotype_list"]) == 1) and (stage_dict["hic_scaffolding"]["parameters"][parameters_label]["option_set"]["use_phased_reads"]):
                    #stage_dict["hic_scaffolding"]["parameters"][parameters_label]["included"] = False
                    stage_dict["hic_scaffolding"]["parameters"].pop(parameters_label)

    #for parameter_label in stage_dict["hic_scaffolding"]["parameters"].keys(): # remove ignore
    #    if not stage_dict["hic_scaffolding"]["parameters"][parameter_label]["included"]:
    #        stage_dict["hic_scaffolding"]["parameters"].pop(parameter_label)

    parameters_list = list(stage_dict["hic_scaffolding"]["parameters"].keys())

    #haplotype=stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"],
    #                            parameters=[parameters_label]) for parameters_label in parameters_list]
    # stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"]
    #prev_parameters_label = stage_dict["hic_scaffolding"]["parameters"][parameters_label]["prev_parameters"]
    if not (config["skip_prescaf_pretext"] or config["skip_both_pretext"]):
        results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{resolution}.map.{ext}",
                                  genome_prefix=[config["genome_prefix"], ],
                                  assembly_stage=[prev_stage,],
                                  haplotype=stage_dict[prev_stage]["parameters"][stage_dict["hic_scaffolding"]["parameters"][current_parameter_label]["prev_parameters"]]["haplotype_list"],
                                  phasing_kmer_length=[stage_dict["hic_scaffolding"]["parameters"][current_parameter_label]["option_set"]["phasing_kmer_length"]], #[stage_dict["hic_scaffolding"]["parameters"][parameters_label]["option_set"]["phasing_kmer_length"] for parameter_label in stage_dict["hic_scaffolding"]["parameters"]],
                                  parameters=[stage_dict["hic_scaffolding"]["parameters"][current_parameter_label]["prev_parameters"]],
                                  resolution=parameters["tool_options"]["pretextsnapshot"]["resolution"],
                                  ext=parameters["tool_options"]["pretextsnapshot"]["format"])  for current_parameter_label in stage_dict["hic_scaffolding"]["parameters"]],
                        ]
    if not (config["skip_postscaf_pretext"] or config["skip_both_pretext"]):
        results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{resolution}.map.{ext}",
                                  genome_prefix=[config["genome_prefix"], ],
                                  assembly_stage=["hic_scaffolding",],
                                  haplotype=stage_dict["hic_scaffolding"]["parameters"][parameters_label]["haplotype_list"],
                                  phasing_kmer_length=[stage_dict["hic_scaffolding"]["parameters"][parameters_label]["option_set"]["phasing_kmer_length"]], #[stage_dict["hic_scaffolding"]["parameters"][parameters_label]["option_set"]["phasing_kmer_length"]],
                                  parameters=[parameters_label],
                                  resolution=parameters["tool_options"]["pretextsnapshot"]["resolution"],
                                  ext=parameters["tool_options"]["pretextsnapshot"]["format"]) for parameters_label in stage_dict["hic_scaffolding"]["parameters"]],
                        ]

    results_list += [
                    [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam.general_stats",
                            genome_prefix=[config["genome_prefix"], ],
                            assembly_stage=[prev_stage,],
                            haplotype=stage_dict[prev_stage]["parameters"][stage_dict["hic_scaffolding"]["parameters"][current_parameter_label]["prev_parameters"]]["haplotype_list"],
                            phasing_kmer_length=[stage_dict["hic_scaffolding"]["parameters"][current_parameter_label]["option_set"]["phasing_kmer_length"]], #[stage_dict["hic_scaffolding"]["parameters"][parameters_label]["option_set"]["phasing_kmer_length"] for parameter_label in stage_dict["hic_scaffolding"]["parameters"]],
                            parameters=[stage_dict["hic_scaffolding"]["parameters"][current_parameter_label]["prev_parameters"]],) for current_parameter_label in stage_dict["hic_scaffolding"]["parameters"]],
                    *[expand(out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["hic_scaffolding", ],
                           haplotype=stage_dict["hic_scaffolding"]["parameters"][parameters_label]["haplotype_list"],
                           parameters=[parameters_label]) for parameters_label in stage_dict["hic_scaffolding"]["parameters"]],
                    [expand(out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.hic",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["hic_scaffolding", ],
                           haplotype=stage_dict["hic_scaffolding"]["parameters"][parameters_label]["haplotype_list"],
                           parameters=[parameters_label]) for parameters_label in stage_dict["hic_scaffolding"]["parameters"]] if not config["skip_hic_file"] else [],
                    expand(out_dir_path / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["hic_scaffolding"],),
                    ]

    for parameters_label in parameters_list:
        if stage_dict["hic_scaffolding"]["parameters"][parameters_label]["hic_scaffolder"] == "yahs":
            if not config["skip_hic_file"]:
                results_list += [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.hic",
                                        genome_prefix=[config["genome_prefix"], ],
                                        assembly_stage=["hic_scaffolding", ],
                                        haplotype=stage_dict["hic_scaffolding"]["parameters"][parameters_label]["haplotype_list"],
                                        parameters=[parameters_label])
                                 ]
    #for parameters_label in parameters_list:
    #    if stage_dict["hic_scaffolding"]["parameters"][parameters_label]["hic_scaffolder"] == "threeddna":
    #        results_list += [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.merged_nodups.txt",
    #                                genome_prefix=[config["genome_prefix"], ],
    #                                assembly_stage=["hic_scaffolding", ],
    #                                haplotype=stage_dict["hic_scaffolding"]["parameters"][parameters_label]["haplotype_list"],
    #                                parameters=[parameters_label])
    #                        ]


    #out_dir_path / "hic_scaffolding/{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.merged_nodups.txt",


    if not config["skip_busco"]:
        results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.busco5.{busco_lineage}.tar.gz",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["hic_scaffolding", ],
                                haplotype=stage_dict["hic_scaffolding"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["hic_scaffolding"]["parameters"]],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["hic_scaffolding"],
                                #haplotype=stage_dict["hic_scaffolding"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["hic_scaffolding"]["parameters"]],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["hic_scaffolding"],
                                haplotype=stage_dict["hic_scaffolding"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["hic_scaffolding"]["parameters"]],
                         expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["hic_scaffolding"],
                                parameters=parameters_list
                                ),
                         ]
#print(stage_dict)
"""

if "gap_closing" in config["stage_list"]: # TODO: modify it and all initiation of stage_dict entries to make it normal!!!!
    stage = "gap_closing"
    prev_stage = stage_dict[stage]["prev_stage"]
    tool_list = config["stage_coretools"][stage]["default"]
    stage_dict[stage]["parameters"] = {}
    for tool in tool_list:
        for option_set in config["coretool_option_sets"][tool]:
            for prev_parameters in ["default"]:
                parameters_label = "{0}_{1}".format(tool, option_set)
                stage_dict[stage]["parameters"][parameters_label] = {}
                stage_dict[stage]["parameters"][parameters_label]["included"] = True
                stage_dict[stage]["parameters"][parameters_label]["gap_closer"] = tool
                stage_dict[stage]["parameters"][parameters_label]["prev_stage"] = prev_stage
                stage_dict[stage]["parameters"][parameters_label]["prev_parameters"] = prev_parameters
                stage_dict[stage]["parameters"][parameters_label]["option_set"] = parameters["tool_options"][tool][option_set] if tool in parameters["tool_options"] else None
                stage_dict[stage]["parameters"][parameters_label]["haplotype_list"] = haplotype_list
"""
if "curation" in config["stage_list"]:
    prev_stage = stage_dict["curation"]["prev_stage"]
    curation_tool_list = config["stage_coretools"]["curation"]["default"]
    stage_dict["curation"]["parameters"] = {}

    for curation_tool in curation_tool_list:
        for option_set in config["coretool_option_sets"][curation_tool]:
            print(prev_stage)
            for prev_parameters in stage_dict[prev_stage]["parameters"]:
                parameters_label = "{0}..{1}_{2}".format(prev_parameters, curation_tool, option_set)
                stage_dict["curation"]["parameters"][parameters_label] = {}
                stage_dict["curation"]["parameters"][parameters_label]["included"] = True
                stage_dict["curation"]["parameters"][parameters_label]["curationeer"] = curation_tool
                stage_dict["curation"]["parameters"][parameters_label]["prev_stage"] = prev_stage
                stage_dict["curation"]["parameters"][parameters_label]["prev_parameters"] = prev_parameters
                stage_dict["curation"]["parameters"][parameters_label]["option_set"] = parameters["tool_options"][curation_tool][option_set] if curation_tool in parameters["tool_options"] else None
                stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"] = stage_dict[stage_dict["curation"]["prev_stage"]]["parameters"][prev_parameters]["haplotype_list"]

    parameters_list = list(stage_dict["curation"]["parameters"].keys())

    if "scaffolds" in  config["curation_seq_type"]:
        if not config["skip_wga"]:
            results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/{target_haplotype}/scaffolds/{genome_prefix}.input.wga.{query_haplotype}.to.{target_haplotype}.YASS.R11.soft.min_len{min_target_len}.png",
                                    genome_prefix=[config["genome_prefix"], ],
                                    assembly_stage=["curation", ],
                                    parameters=[parameters_label],
                                    min_target_len=parameters["tool_options"]["wga"][ "min_target_len"],
                                    query_haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"] + list(input_reference_filedict.keys()),
                                    target_haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                    ) for parameters_label in stage_dict["curation"]["parameters"]]]
        results_list += [[[[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/scaffolds/{genome_prefix}.input.{haplotype}.{track_type}.win{window}.step{step}.{threshold_type}.png",
                                threshold_type=["absolute", "relative"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                track_type=[track_type],
                                window=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["window"]],
                                step=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["step"]],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for window_settings in stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"] ] for parameters_label in stage_dict["curation"]["parameters"]] for track_type in ("gap", "windowmasker", "trf", "gc") ],
                         [[[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{track_type}.win{window}.step{step}.track.stat",
                                seq_type=["scaffolds"],
                                threshold_type=["absolute", "relative"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                track_type=[track_type],
                                window=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["window"]],
                                step=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["step"]],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for window_settings in stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"] ] for parameters_label in stage_dict["curation"]["parameters"]] for track_type in ("gap", "windowmasker", "trf", "gc") ],
                         [[expand(out_dir_path / "curation/{parameters}/{haplotype}/scaffolds/{genome_prefix}.input.{haplotype}.{datatype}.coverage.win{window}.step{step}.png",
                                window=stage_dict["curation"]["parameters"][parameters_label]["option_set"]["coverage"]["options"][window_step_set]["window"],
                                step=stage_dict["curation"]["parameters"][parameters_label]["option_set"]["coverage"]["options"][window_step_set]["step"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                datatype=coverage_track_data_type_set,
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for window_step_set in stage_dict["curation"]["parameters"][parameters_label]["option_set"]["coverage"]["options"]] for parameters_label in stage_dict["curation"]["parameters"]] if coverage_track_data_type_set else [],
                         [[expand(out_dir_path / "curation/{parameters}/{haplotype}/scaffolds/{genome_prefix}.input.{haplotype}.{datatype}_{cov_type}_coverage.win{window}.step{step}.track.bedgraph",
                                cov_type=["mean", "median"],
                                window=stage_dict["curation"]["parameters"][parameters_label]["option_set"]["coverage"]["options"][window_step_set]["window"],
                                step=stage_dict["curation"]["parameters"][parameters_label]["option_set"]["coverage"]["options"][window_step_set]["step"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                datatype=coverage_track_data_type_set,
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for window_step_set in stage_dict["curation"]["parameters"][parameters_label]["option_set"]["coverage"]["options"]] for parameters_label in stage_dict["curation"]["parameters"]] if coverage_track_data_type_set else [],
                         [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.canonical.txt",
                                seq_type=["scaffolds"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],
                         [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.cannonical_telomere.win1000.step200.track.bedgraph",
                                seq_type=["scaffolds"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],
                         [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.cannonical_telomere_warning.win1000.step200.track.bedgraph",
                                seq_type=["scaffolds"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],
                         [expand(out_dir_path  / "{assembly_stage}/{parameters}/{haplotype}/scaffolds/{genome_prefix}.input.{haplotype}.{datatype}.vcf.gz",
                                genome_prefix=[config["genome_prefix"], ],
                                datatype=variant_calling_data_type_set,
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]] if variant_calling_data_type_set and (not config["skip_variantcalling"]) else [],
                         ]
        if "hic_scaffolding" in config["stage_list"]:
            results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.cannonical_telomere_warning.win1000.step200.scaled.track.bedgraph",
                                seq_type=["scaffolds"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],
                            [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.cannonical_telomere.win1000.step200.scaled.track.bedgraph",
                                seq_type=["scaffolds"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],
                            [[expand(out_dir_path / "curation/{parameters}/{haplotype}/scaffolds/{genome_prefix}.input.{haplotype}.{datatype}_{cov_type}_coverage.win{window}.step{step}.scaled.track.bedgraph",
                                cov_type=["mean", "median"],
                                window=stage_dict["curation"]["parameters"][parameters_label]["option_set"]["coverage"]["options"][window_step_set]["window"],
                                step=stage_dict["curation"]["parameters"][parameters_label]["option_set"]["coverage"]["options"][window_step_set]["step"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                datatype=coverage_track_data_type_set,
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for window_step_set in stage_dict["curation"]["parameters"][parameters_label]["option_set"]["coverage"]["options"]] for parameters_label in stage_dict["curation"]["parameters"]] if coverage_track_data_type_set else [],
                            [[[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{track_type}.win{window}.step{step}.scaled.track.bedgraph",
                                seq_type=["scaffolds"],
                                threshold_type=["absolute", "relative"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                track_type=[track_type],
                                window=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["window"]],
                                step=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["step"]],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for window_settings in stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"] ] for parameters_label in stage_dict["curation"]["parameters"]] for track_type in ("gap", "windowmasker", "trf", "gc") ],
                             ]

    if ("contigs" in config["curation_seq_type"]) and (prev_stage == "hic_scaffolding") :
        results_list += [[[[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{track_type}.win{window}.step{step}.track.{filetype}",
                                filetype=["stat", "assembly.bedgraph"],
                                seq_type=["contigs"],
                                threshold_type=["absolute", "relative"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                track_type=[track_type],
                                window=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["window"]],
                                step=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["step"]],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for window_settings in stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"] ] for parameters_label in stage_dict["curation"]["parameters"]] for track_type in ("gap", "windowmasker", "trf", "gc") ],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.canonical.txt",
                                seq_type=["contigs"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{telo_type}.win1000.step200.track.assembly.bedgraph",
                                telo_type=["cannonical_telomere", "non_cannonical_telomere"],
                                seq_type=["contigs"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{telo_type}.win1000.step200.track.assembly.bedgraph",
                                telo_type=["cannonical_telomere_warning", "non_cannonical_telomere_warning"],
                                seq_type=["contigs"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],
                         ]

    if not config["skip_higlass"]:
        for parameters_label in parameters_list:
            if stage_dict["curation"]["parameters"][parameters_label]["prev_stage"] == "hic_scaffolding": # TODO: add handling for a case when "hic_scaffolding" is not a stage before the "curation"
                results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/scaffolds/{genome_prefix}.input.{haplotype}.higlass.mcool",
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],]


    if config["create_hic_file_during_curation"]:
        results_list += [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.pre.hic",
                                assembly_stage=[stage_dict["curation"]["prev_stage"]],
                                haplotype=stage_dict[prev_stage]["parameters"][stage_dict["curation"]["parameters"][current_parameter_label]["prev_parameters"]]["haplotype_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                phasing_kmer_length=["31"],# [stage_dict["curation"]["parameters"][current_parameter_label]["option_set"]["phasing_kmer_length"]],
                                parameters=[stage_dict["curation"]["parameters"][current_parameter_label]["prev_parameters"]],
                                ) for current_parameter_label in stage_dict["curation"]["parameters"]]

#----

#---- Final rule ----
rule all:
    input:
        results_list
        #results_dict[config["mode"]]
#----

#---- Include section ----
include: "workflow/rules/Install/Pip.smk"
include: "workflow/rules/Preprocessing/Files.smk"
include: "workflow/rules/QCFiltering/FastQC.smk"
include: "workflow/rules/QCFiltering/MultiQC.smk"
include: "workflow/rules/QCFiltering/Cutadapt.smk"
include: "workflow/rules/QCFiltering/Trimmomatic.smk"

#if "nanopore" in data_types:
include: "workflow/rules/QCFiltering/Nanopore.smk"
include: "workflow/rules/QCFiltering/NanoQC.smk"
include: "workflow/rules/QCFiltering/NanoPlot.smk"

include: "workflow/rules/Kmer/Jellyfish.smk"
include: "workflow/rules/Kmer/Meryl.smk"
include: "workflow/rules/Kmer/Smudgeplot.smk"
#include: "workflow/rules/Kmer/KAT.smk"
include: "workflow/rules/Kmer/GCplot.smk"
include: "workflow/rules/Kmer/Genomescope.smk"

include: "workflow/rules/QCAssembly/BUSCO5.smk"
include: "workflow/rules/QCAssembly/Merqury.smk"
include: "workflow/rules/QCAssembly/QUAST.smk"
include: "workflow/rules/QCAssembly/General.smk"
include: "workflow/rules/Contamination/FCS.smk"
include: "workflow/rules/Contamination/Kraken2.smk"

if "hifi" in data_types:
    include: "workflow/rules/Contigs/Hifiasm.smk"

include: "workflow/rules/Contigs/Graph.smk"
include: "workflow/rules/Stats/General.smk"

if "purge_dups" in config["stage_list"]:
    include: "workflow/rules/Purge_dups/Purge_dups.smk"
    include: "workflow/rules/Purge_dups/Purge_dupsQC.smk"

include: "workflow/rules/HiC/ReadPhasing.smk"

include: "workflow/rules/Alignment/Index.smk"
include: "workflow/rules/Alignment/Stats.smk"

if ("hic_scaffolding" in config["stage_list"]) or ("curation" in config["stage_list"]) or ("gap_closing" in config["stage_list"]):
    if config["other_tool_option_sets"]["mapping_pipeline"] == "arima":
        include: "workflow/rules/Alignment/Arima.smk"
    elif config["other_tool_option_sets"]["mapping_pipeline"] == "bwa_only":
        include: "workflow/rules/Alignment/BWAOnly.smk"
    elif config["other_tool_option_sets"]["mapping_pipeline"] == "pairtools":
        include: "workflow/rules/Alignment/Pairtools.smk"
    include: "workflow/rules/Alignment/PostAlignment.smk"

if ("hic_scaffolding" in config["stage_list"]) or ("curation" in config["stage_list"]):
    include: "workflow/rules/Alignment/Pretext.smk"

if "hic_scaffolding" in config["stage_list"]:
    #include: "workflow/rules/HiC/Salsa2.smk"
    include: "workflow/rules/HiC/YAHS.smk"
    include: "workflow/rules/HiC/3DDNA.smk"

if "curation" in config["stage_list"]:
    include: "workflow/rules/Curation/RapidCuration.smk"
    include: "workflow/rules/Curation/GapTrack.smk"
    include: "workflow/rules/Curation/WindowmaskerTrack.smk"
    include: "workflow/rules/Curation/CoverageTrack.smk"
    include: "workflow/rules/Curation/TelomereTrack.smk"
    include: "workflow/rules/Curation/HiGlassTrack.smk"
    include: "workflow/rules/Curation/TRFTrack.smk"
    include: "workflow/rules/Curation/Masking.smk"
    include: "workflow/rules/Curation/GCTrack.smk"
    include: "workflow/rules/Curation/WGA.smk"
    include: "workflow/rules/Curation/VariantTrack.smk"

if "gap_closing" in config["stage_list"]:
    include: "workflow/rules/Finalization/GapClosing.smk"