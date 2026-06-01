import os
import sys
from functools import partial
import yaml
#import logging
import shutil
from copy import deepcopy
from collections import OrderedDict
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path, PosixPath
from numbers import Number # Abstract class for numeric types
import pandas as pd
from stone.backends.obj_c import comment_prefix

#---- Include sections for functions ----
include: "workflow/functions/option_parsing.py"
include: "workflow/functions/general_parsing.py"
include: "workflow/functions/resources.py"
#----------------------------------------

#---- Read config files ----
#-------- Read core config file --------
with open(config["main_config_file"], "r") as core_yaml_fd:
    config.update(yaml.safe_load(core_yaml_fd))
#-------- Read secondary tools config file -------
with open(config["secondary_tool_config_file"], "r") as secondary_tool_fd:
    copy_absent_entries(yaml.safe_load(secondary_tool_fd), config)
#-------- Read cluster config file ---------
with open(config["cluster_config_file"], "r") as cluster_fd:
    copy_absent_entries(yaml.safe_load(cluster_fd), config)
#-------- Read 'skip' config file --------
with open(config["skip_config_file"], "r") as skip_yaml_fd:
    for key, value in yaml.safe_load(skip_yaml_fd).items():
        if key not in config:
            config[key] = value
#---------------------------------------
#-------- Read resources config files --------
resources_dir_path = Path(config["resources_dir"])
for resource, res_datatype in zip(["threads", "memory_mb", "time"], [int, int, str]):
    resource_df = pd.read_csv(resources_dir_path / f"{config['resource_profile']}/{resource}.tab",
                              sep="\t", header=0, index_col=0)
    for config_label in resource_df.columns:
        config["parameters"][config_label][resource] = resource_df[config_label].to_dict(OrderedDict)

#---------------------------------------------

#---------------------------




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
for datatype in config["final_kmer_datatypes"]:
    if datatype not in fastq_based_data_type_set:
        if config["mode"] in ["preprocessing", "qc"]:
            pass
        else:
            if ("skip_kmer" in config) and (config["skip_kmer"]):
                pass
            else:
                raise ValueError("ERROR!!! dinal kmer datatype ({0}) is absent among input fastq-based datatypes({1})".format(datatype,
                                                                                                                              ",".join(fastq_based_data_type_set)))

#--------

#----

#---- Checking input files ----
candidate_agp_dir_path = input_dir_path / "candidate_chr/"

candidate_agp_filename = list(candidate_agp_dir_path.glob("*.agp"))
#print(candidate_agp_filename)
#print(candidate_agp_dir_path.name)
if len(candidate_agp_filename) > 1:
    raise ValueError(f"ERROR!!! More than one agp file was detected in folder {str(candidate_agp_dir_path.name)}!")
elif len(candidate_agp_filename) == 1:
    candidate_agp_filename = candidate_agp_filename[0]
    candidate_output_dir = output_dict["data"] / "candidate_chr/"
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
    #print(chr_component_series)
    candidate_chr_id_list = list(chr_component_series.index)
    for scaffold_id in chr_component_series.index:
        #print(chr_component_series[scaffold_id])
        chr_component_series[[scaffold_id]].to_csv(f"{candidate_output_prefix}.{scaffold_id}.components.ids",sep="\t",header=False,index=False)
        chr_black_list_series = chr_component_series[~chr_component_series.isin(chr_component_series[[scaffold_id]])]
        chr_black_list_series.to_csv(f"{candidate_output_prefix}.{scaffold_id}.pretext.blacklist",sep="\t",header=False,index=False)
else:
    candidate_agp_filename = []
#logging.info("Checking input files...")

input_filedict = OrderedDict()
input_file_prefix_dict = OrderedDict()
input_fasta_filedict = OrderedDict()
input_fasta_file_prefix_dict = OrderedDict()
input_reference_filedict = OrderedDict()

input_forward_suffix_dict = OrderedDict()
input_reverse_suffix_dict = OrderedDict()
input_pairprefix_dict = OrderedDict()

for d_type in fastq_based_data_type_set:
    input_filedict[d_type] = find_fastqs(input_dict[d_type]["fastq_dir"], fastq_extension=config["fastq_extension"])
    input_file_prefix_dict[d_type] = list(map(lambda s: str(s.name)[:-len(config["fastq_extension"])],
                                                input_filedict[d_type]))

for d_type in fasta_based_data_type_set:
    input_fasta_filedict[d_type] = find_fastas(input_dict[d_type]["fasta_dir"], fasta_extension=config["fasta_extension"])
    input_fasta_file_prefix_dict[d_type] = list(map(lambda s: str(s.name)[:-len(config["fasta_extension"])],
                                                input_fasta_filedict[d_type]))

#---- detect datatypes and check if datatype has files in both fasta and fastq formats ----
datatype_format_dict = OrderedDict()
datatype_extension_dict = OrderedDict()
for d_type in set(data_types):
    print(d_type)
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
    reference_input_dir = input_dict["reference"]["dir"]
    reference_genomes_list = []
    for element in reference_input_dir.glob("*"):
        if element.is_dir():
            reference_genomes_list.append(element.name)
    print(reference_input_dir)
    print(reference_genomes_list)
    for genome in reference_genomes_list:
        input_reference_filedict[genome] = {}
        for filetype in "fasta", "syn", "whitelist", "orderlist":
            input_reference_filedict[genome][filetype] = list((reference_input_dir / genome).glob("*.{0}".format(filetype)))
            print(input_reference_filedict[genome][filetype])
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

for tool in parameters["tool_options"]: # sort datatypes in case of mixed datatypes to avoid double calculations
    for option_set in parameters["tool_options"][tool]:
        if "main_datatypes" in parameters["tool_options"][tool][option_set]:
            parameters["tool_options"][tool][option_set]["main_datatypes"] = sorted(parameters["tool_options"][tool][option_set]["main_datatypes"])
        if "qc_datatypes" in parameters["tool_options"][tool][option_set]:
            parameters["tool_options"][tool][option_set]["qc_datatypes"] = sorted(parameters["tool_options"][tool][option_set]["qc_datatypes"])


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
#include: "workflow/rules/Contigs/Hifiasm.smk"
print(stage_dict)
print(datatype_format_dict)
print(input_file_prefix_dict)
print(input_fasta_file_prefix_dict)
results_dict = {}

haplotype_list = ["hap{0}".format(i) for i in range(1, config["ploidy"] + 1)] # TODO: obsolete: remove and fix issues
primary_haplotype = "hap1" # TODO: obsolete: remove and fix issues

results_list = []

"""
# Pip dependencies are now installed azutomatically
for conda_env in config["conda"]:
    if "pip" in config["conda"][conda_env]:
        results_list += [expand("results/config/pip.{conda_env}.requirements",
                                conda_env=[conda_env])]
"""

#---- Create output filelist ----
if "check_reads" in config["stage_list"]:
    results_list += [
                     final_config_yaml,
                     final_input_yaml
                     ]

if "check_draft" in config["stage_list"]:
    results_list += [ ] # TODO: implement

if ("read_qc" in config["stage_list"]) and (not config["skip_read_qc"]):
    results_list += [[expand(output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip",
                               datatype=[dat_type, ],
                               stage=["raw", ],
                               fileprefix=input_file_prefix_dict[dat_type],) for dat_type in fastqc_data_type_set ],
                     expand(output_dict["qc"] / "fastqc/hic/{stage}/{fileprefix}_fastqc.zip",
                         stage=["orig", ],
                         fileprefix=input_file_prefix_dict["hic"],) if ("hic" in fastqc_data_type_set) and (config["hic_enzyme_set"] == 'Arima') else [] ,
                     expand(output_dict["qc"] / "multiqc/{datatype}/{stage}/multiqc.{datatype}.{stage}.report.html",
                             datatype=fastqc_data_type_set ,
                             stage=["raw",]),
                     #[expand(output_dict["qc"] / "nanoplot/{datatype}/{stage}/{fileprefix}.Yield_By_Length.png",
                     #          datatype=[dat_type, ],
                     #          stage=["raw", ],
                     #          fileprefix=input_file_prefix_dict[dat_type],) for dat_type in long_read_data_type_set],
                     [expand(output_dict["qc"] / "nanoplot/{datatype}/{stage}/{datatype}.{stage}.NanoStats.tsv",
                               datatype=[dat_type, ],
                               stage=["raw", ],
                               ) for dat_type in long_read_data_type_set],
                     [expand(output_dict["qc"] / "nanoqc/{datatype}/{stage}/{fileprefix}",
                               datatype=[dat_type, ],
                               stage=["raw", ],
                               fileprefix=input_file_prefix_dict[dat_type],) for dat_type in long_read_data_type_set],
                     ]
    if ("hic" in data_types) and ((config["hic_enzyme_set"] == "custom") or config["hic_enzyme_dict"][config["hic_enzyme_set"]]):
        #print(data_types)
        results_list += [expand(output_dict["qc"] / "tadbit/hic/raw/{genome_prefix}.tadbit.stats",
            genome_prefix=[config["genome_prefix"]])]

if not config["skip_mtdna"]:
    if not config["skip_mitohifi_reads"]:
        if not config["skip_mitohifi_reads_per_file"]:
            if "hifi" in data_types:
                results_list += [ expand(out_dir_path / "mtDNA/mitohifi/{mtdna_ref}/hifi/filtered/{fileprefix}/FINISH_FLAG",
                                         mtdna_ref=["recommended"],
                                         fileprefix=input_file_prefix_dict["hifi"])]
        if not config["skip_mitohifi_reads_combined"]:
            if "hifi" in data_types:
                results_list += [expand(out_dir_path / "mtDNA/mitohifi/{mtdna_ref}/hifi/combined/hifi.combined/FINISH_FLAG",
                                        mtdna_ref=["recommended"],)]
    if not config["skip_mitoz"]:
        if (not config["skip_mitoz_hic"]) and ("hic" in data_types):
            results_list += [expand(out_dir_path / "mtDNA/mitoz/denovo/{datatype}/{stage}/{pairprefix}/FINISH_FLAG",
                                    datatype=["hic",],
                                    stage=["filtered"],
                                    pairprefix=input_pairprefix_dict["hic"])]
        if (not config["skip_mitoz_illumina"]) and ("illumina" in data_types):
            results_list += [expand(out_dir_path / "mtDNA/mitoz/denovo/{datatype}/{stage}/{pairprefix}/FINISH_FLAG",
                                    datatype=["illumina",],
                                    stage=["filtered"],
                                    pairprefix=input_pairprefix_dict["illumina"])]

if "draft_qc" in config["stage_list"]:
    current_stage = "draft_qc"
    draft_file_dict = get_input_assemblies(input_dir_path / "draft/fasta", config["ploidy"], config["assembly_fasta_extension"])
    #print(draft_file_dict)
    stage_dict["draft_qc"]["parameters"] = {}

    for qcer in config["stage_coretools"]["draft_qc"]["default"]:
        for option_set in config["coretool_option_sets"][qcer]:
            parameters_label="{0}_{1}".format(qcer, option_set)
            stage_dict["draft_qc"]["parameters"][parameters_label] = {}
            stage_dict["draft_qc"]["parameters"][parameters_label]["qcer"] = qcer
            stage_dict["draft_qc"]["parameters"][parameters_label]["stage_seq_type"] = None
            #print(qcer)
            stage_dict["draft_qc"]["parameters"][parameters_label]["option_set"] = deepcopy(parameters["tool_options"][qcer]) #[option_set])
            stage_dict["draft_qc"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] = config["ploidy"]
            stage_dict["draft_qc"]["parameters"][parameters_label]["haplotype_list"] = ["hap{0}".format(i) for i in range(1, stage_dict["draft_qc"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] + 1)] if stage_dict["draft_qc"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] > 1 else ["hap0"]
            stage_dict["draft_qc"]["parameters"][parameters_label]["option_set_group"] = None

            #TEMP STRING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO: replace
            stage_dict["draft_qc"]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"] = ["hifi"]

            #if not stage_dict["draft_qc"]["parameters"][parameters_label]["option_set"]["qc_datatypes"]:
            #    pass # TODO: ADD!!!

            #if ("purge_dups" in config["qc_settings"]) and ("qc_datatypes" in config["qc_settings"]["purge_dups"]) and config["qc_settings"]["purge_dups"]["qc_datatypes"]:
            #    stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"] = config["qc_settings"]["purge_dups"]["qc_datatypes"]
            #else:
            #
            #    #if not stage_dict["contig"]["parameters"][parameters_label]["option_set"]["qc_datatypes"]:
            #    #    stage_dict["contig"]["parameters"][parameters_label]["option_set"]["qc_datatypes"] = stage_dict["contig"]["parameters"][parameters_label]["option_set"]["main_datatypes"]
            #    if not stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"]:
            #        stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"] = stage_dict["contig"]["parameters"][parameters_label]["option_set"]["main_datatypes"]



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
    results_list += [
                     [[[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/trackplots/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.{scaffold_length}.win{window}.step{step}.{threshold_type}.png",
                               scaffold_length=config["qc_settings"]["assembly_scaffold_sets"],
                               threshold_type=config["qc_settings"]["threshold_types"],
                               genome_prefix=[config["genome_prefix"], ],
                               assembly_stage=[current_stage, ],
                               track_type=[track_type],
                               window=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["window"]],
                               step=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["step"]],
                               haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                               parameters=[parameters_label])

                        for window_settings in config["qc_settings"]["windows_sets"]]
                        for parameters_label in stage_dict[current_stage]["parameters"]]
                        for track_type in ["gc"]],  #"windowmasker", "trf"
                     [expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.bedgraph",
                            genome_prefix=[config["genome_prefix"], ],
                            assembly_stage=[current_stage, ],
                            parameters=[parameters_label],
                            haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                           ) for parameters_label in stage_dict[current_stage]["parameters"]
                     ]
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

    if (not config["skip_all_pretext"]) and (not config["skip_draft_qc_pretext"]):
        results_list += [[[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.NA.{haplotype}.{subset}.rmdup.mapq{mapq}.{res}.tracks.pretext",
                                  res=["high_res"], #"default",
                                  haplotype=["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                                  subset=["all"] , # + (["microchr"] if ("bird_genome" in config) and config["bird_genome"] else [])
                                  genome_prefix=[config["genome_prefix"], ],
                                  assembly_stage=[current_stage],
                                  parameters=stage_dict[current_stage]["parameters"],
                                  resolution=parameters["tool_options"]["pretextsnapshot"]["resolution"],
                                  mapq=parameters["tool_options"]["pretextmap"]["mapq"],
                                  ext=parameters["tool_options"]["pretextsnapshot"]["format"],
                                  #window=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["window"],
                                  #step = parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["step"],
                                     ) for window_step_set in config["qc_settings"]["windows_sets"]] for parameters_label in stage_dict[current_stage]["parameters"]] if coverage_track_data_type_set else [],
                         ]
    if (not config["skip_wga"]) and (not config["skip_draft_qc_wga"]):
        results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/wga.{query_prefix}.{query_length}.to.{target_prefix}.{target_length}.YASS.R11.soft.min_len{min_target_len}.png",
                                     query_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     target_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     genome_prefix=[config["genome_prefix"], ],
                                     assembly_stage=[current_stage, ],
                                     parameters=[parameters_label],
                                     min_target_len=parameters["tool_options"]["wga"]["min_target_len"],
                                     query_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                     target_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                   ) for parameters_label in stage_dict[current_stage]["parameters"]],
                             ]
        if input_reference_filedict:
            results_list += [
                             [expand(out_dir_path / "{assembly_stage}/{parameters}/wga.{query_prefix}.{query_length}.to.{target_prefix}.{target_length}.YASS.R11.soft.min_len{min_target_len}.png",
                                     query_length=config["qc_settings"]["reference_scaffold_sets"],
                                     target_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     genome_prefix=[config["genome_prefix"], ],
                                     assembly_stage=[current_stage, ],
                                     parameters=[parameters_label],
                                     min_target_len=parameters["tool_options"]["wga"]["min_target_len"],
                                     query_prefix=list(input_reference_filedict.keys()),
                                     target_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                   ) for parameters_label in stage_dict[current_stage]["parameters"]],
                             ]

    if input_reference_filedict and (not config["skip_draft_qc_ragtag"]) and (not config["skip_ragtag"]):
        results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/ragtag/{reference}/{genome_prefix}.{assembly_stage}.{haplotype}.to.{reference}.fasta",
                                        genome_prefix=[config["genome_prefix"], ],
                                        assembly_stage=[current_stage],
                                        parameters=[parameters_label],
                                        reference=list(input_reference_filedict.keys()),
                                        haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                                        ) for parameters_label in stage_dict[current_stage]["parameters"]],
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

    if candidate_agp_filename:
        results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/per_chr/{genome_prefix}.{assembly_stage}.{haplotype}.NA.{candidate_chr_id}.rmdup.precurated.mapq{mapq}.{res}.tracks.pretext",
                                candidate_chr_id=candidate_chr_id_list,
                                assembly_stage=[current_stage],
                                parameters=[parameter_label],
                                haplotype=["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                                genome_prefix=[config["genome_prefix"], ],
                                res=["high_res"],
                                mapq=parameters["tool_options"]["pretextmap"]["mapq"],
                                #window=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["window"],
                                #step=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["step"],
                                ) for parameter_label in stage_dict[current_stage]["parameters"] for window_step_set in config["qc_settings"]["windows_sets"]]
                         ]

    if "gap_closing" in config["stage_list"]:
        current_stage = "gap_closing"
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
                    if not stage_dict["gap_closing"]["parameters"][parameters_label]["option_set"]["qc_datatypes"]:
                        stage_dict["gap_closing"]["parameters"][parameters_label]["option_set"]["qc_datatypes"] = stage_dict["gap_closing"]["parameters"][parameters_label]["option_set"]["main_datatypes"]
                    #TEMP STRING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TODO: replace
                    stage_dict["gap_closing"]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"] = ["hifi"]

        parameters_list = list(stage_dict["gap_closing"]["parameters"].keys())
        results_list += [*[expand(out_dir_path / "gap_closing/{parameters}/{genome_prefix}.gap_closing.{haplotype}.len",
                                parameters=[parameters_label],
                                genome_prefix=[config["genome_prefix"], ],
                                haplotype=stage_dict["gap_closing"]["parameters"][parameters_label]["haplotype_list"]
                                ) for parameters_label in parameters_list],
                         expand(out_dir_path / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["gap_closing"],),
                        ]
        if not config["skip_hic_file"]:
            results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.pre.mapq{mapq}.hic",
                                      assembly_stage=["gap_closing"],
                                      parameters=[parameters_label],
                                      genome_prefix=[config["genome_prefix"], ],
                                      haplotype=stage_dict["gap_closing"]["parameters"][parameters_label]["haplotype_list"],
                                      phasing_kmer_length=[stage_dict["gap_closing"]["parameters"][parameters_label]["option_set"]["phasing_kmer_length"]],
                                      mapq=parameters["tool_options"]["yahs_juicer_pre"]["mapq"])
                           for parameters_label in parameters_list] if not config["skip_hic_file"] else []
                         ]

        results_list += [
                     [[[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/trackplots/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.{scaffold_length}.win{window}.step{step}.{threshold_type}.png",
                               scaffold_length=config["qc_settings"]["assembly_scaffold_sets"],
                               threshold_type=config["qc_settings"]["threshold_types"],
                               genome_prefix=[config["genome_prefix"], ],
                               assembly_stage=[current_stage, ],
                               track_type=[track_type],
                               window=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["window"]],
                               step=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["step"]],
                               haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                               parameters=[parameters_label])

                        for window_settings in config["qc_settings"]["windows_sets"]]
                        for parameters_label in stage_dict[current_stage]["parameters"]]
                        for track_type in  ["gc"]],  #"windowmasker", "trf"
                     [expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.bedgraph",
                            genome_prefix=[config["genome_prefix"], ],
                            assembly_stage=[current_stage, ],
                            parameters=[parameters_label],
                            haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                           ) for parameters_label in stage_dict[current_stage]["parameters"]
                     ]
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
        if (not config["skip_all_pretext"]) and (not config["skip_gap_closing_pretext"]):
            results_list += [[[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.NA.{haplotype}.{subset}.rmdup.mapq{mapq}.{res}.tracks.pretext",
                                  res=["high_res"], #"default",
                                  haplotype=["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                                  subset=["all"] , # + (["microchr"] if ("bird_genome" in config) and config["bird_genome"] else [])
                                  genome_prefix=[config["genome_prefix"], ],
                                  assembly_stage=[current_stage],
                                  parameters=stage_dict[current_stage]["parameters"],
                                  resolution=parameters["tool_options"]["pretextsnapshot"]["resolution"],
                                  mapq=parameters["tool_options"]["pretextmap"]["mapq"],
                                  ext=parameters["tool_options"]["pretextsnapshot"]["format"],
                                  #window=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["window"],
                                  #step = parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["step"],
                                 ) for window_step_set in config["qc_settings"]["windows_sets"]] for parameters_label in
                                                          stage_dict[current_stage]["parameters"]] if coverage_track_data_type_set else [],
                     ]
        if candidate_agp_filename:
            results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/per_chr/{genome_prefix}.{assembly_stage}.{haplotype}.NA.{candidate_chr_id}.rmdup.precurated.mapq{mapq}.{res}.tracks.pretext",
                                    candidate_chr_id=candidate_chr_id_list,
                                    assembly_stage=[current_stage],
                                    parameters=[parameter_label],
                                    haplotype=["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                                    genome_prefix=[config["genome_prefix"], ],
                                    res=["high_res"],
                                    mapq=parameters["tool_options"]["pretextmap"]["mapq"],
                                    #window=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["window"],
                                    #step=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["step"],
                                    ) for parameter_label in stage_dict[current_stage]["parameters"] for window_step_set in config["qc_settings"]["windows_sets"]]
                             ]

        if (not config["skip_wga"]) and (not config["skip_gap_closing_wga"]):
            results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/wga.{query_prefix}.{query_length}.to.{target_prefix}.{target_length}.YASS.R11.soft.min_len{min_target_len}.png",
                                         query_length=config["qc_settings"]["assembly_scaffold_sets"],
                                         target_length=config["qc_settings"]["assembly_scaffold_sets"],
                                         genome_prefix=[config["genome_prefix"], ],
                                         assembly_stage=[current_stage, ],
                                         parameters=[parameters_label],
                                         min_target_len=parameters["tool_options"]["wga"]["min_target_len"],
                                         query_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                             genome_prefix=[config["genome_prefix"], ],
                                                             assembly_stage=[current_stage, ],
                                                             haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                         target_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                             genome_prefix=[config["genome_prefix"], ],
                                                             assembly_stage=[current_stage, ],
                                                             haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                       ) for parameters_label in stage_dict[current_stage]["parameters"]],
                                 ]
            if input_reference_filedict:
                results_list += [
                                 [expand(out_dir_path / "{assembly_stage}/{parameters}/wga.{query_prefix}.{query_length}.to.{target_prefix}.{target_length}.YASS.R11.soft.min_len{min_target_len}.png",
                                         query_length=config["qc_settings"]["reference_scaffold_sets"],
                                         target_length=config["qc_settings"]["assembly_scaffold_sets"],
                                         genome_prefix=[config["genome_prefix"], ],
                                         assembly_stage=[current_stage, ],
                                         parameters=[parameters_label],
                                         min_target_len=parameters["tool_options"]["wga"]["min_target_len"],
                                         query_prefix=list(input_reference_filedict.keys()),
                                         target_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                             genome_prefix=[config["genome_prefix"], ],
                                                             assembly_stage=[current_stage, ],
                                                             haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                       ) for parameters_label in stage_dict[current_stage]["parameters"]],
                                 ]

        if input_reference_filedict and (not config["skip_gap_closing_ragtag"]) and (not config["skip_ragtag"]):
            results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/ragtag/{reference}/{genome_prefix}.{assembly_stage}.{haplotype}.to.{reference}.fasta",
                                            genome_prefix=[config["genome_prefix"], ],
                                            assembly_stage=[current_stage],
                                            parameters=[parameters_label],
                                            reference=list(input_reference_filedict.keys()),
                                            haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                                            ) for parameters_label in stage_dict[current_stage]["parameters"]],]


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

                    [[expand(output_dict["kmer"] / "{datatype}/{stage}/{analysis_tool}/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{analysis_tool}.parameters",
                           datatype=[dat_type,],
                           genome_prefix=[config["genome_prefix"], ],
                           analysis_tool=["genomescope"] + (["krater"] if not config["skip_krater"] else []),
                           stage=["filtered",],
                           kmer_tool=[kmer_tool,],
                           kmer_length=parameters["tool_options"][kmer_tool][dat_type]["kmer_length"],
                           ) for kmer_tool in config["kmer_counter_list"] ]  for dat_type in genome_size_estimation_data_type_set],


                    ]
    if not config["skip_per_lib_genome_estimation"]:
        results_list += [[expand(output_dict["kmer"] / "{datatype}/{stage}/{analysis_tool}/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{read_prefix}.{analysis_tool}.parameters",
                           datatype=[dat_type,],
                           genome_prefix=[config["genome_prefix"], ],
                           analysis_tool=["genomescope"] + (["krater"] if not config["skip_krater"] else []),
                           stage=["filtered",],
                           kmer_tool=[kmer_tool,],
                           kmer_length=parameters["tool_options"][kmer_tool][dat_type]["kmer_length"],
                           read_prefix=input_pairprefix_dict[dat_type] if dat_type in config["paired_fastq_based_data"] else input_file_prefix_dict[dat_type],
                           ) for kmer_tool in config["kmer_counter_list"] ]  for dat_type in genome_size_estimation_data_type_set],

    results_list += [expand(output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip",
                            datatype=[dat_type, ],
                            stage=["filtered", ],
                            fileprefix=[pairprefix + suffix for suffix in ("_1", "_2") for pairprefix in input_pairprefix_dict[dat_type]],
                            ) for dat_type in set(config["paired_fastq_based_data"]) & fastq_based_data_type_set ]

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
                        expand(output_dict["qc"] / "nanoqc/{datatype}/{stage}/{fileprefix}",
                            datatype=["simplex", ],
                            stage=["trimmed", ],
                            fileprefix=input_file_prefix_dict["simplex"],) if "simplex" in long_read_data_type_set else [],
                        expand(output_dict["qc"] / "nanoqc/{datatype}/{stage}/{fileprefix}",
                            datatype=["duplex", ],
                            stage=["trimmed", ],
                            fileprefix=input_file_prefix_dict["duplex"],) if "duplex" in long_read_data_type_set else [],

                        ]
    if not config["skip_nanoplot"]:
        results_list += [[expand(output_dict["qc"] / "nanoplot/{datatype}/{stage}/{datatype}.{stage}.NanoStats.tsv",
                               datatype=[dat_type, ],
                               stage=["filtered", ],
                               ) for dat_type in long_read_data_type_set],
                        expand(output_dict["qc"] / "nanoplot/{datatype}/{stage}/{datatype}.{stage}.NanoStats.tsv",
                               datatype=["nanopore", ],
                               stage=["trimmed", ],
                               ) if "nanopore" in long_read_data_type_set else [],
                        expand(output_dict["qc"] / "nanoplot/{datatype}/{stage}/{datatype}.{stage}.NanoStats.tsv",
                             datatype=["duplex", ],
                             stage=["trimmed", ],
                                ) if "duplex" in long_read_data_type_set else [],
                        expand(output_dict["qc"] / "nanoplot/{datatype}/{stage}/{datatype}.{stage}.NanoStats.tsv",
                             datatype=["simplex", ],
                             stage=["trimmed", ],
                                ) if "simplex" in long_read_data_type_set else [],
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
    current_stage = "contig"
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
            stage_dict["contig"]["parameters"][parameters_label]["option_set_group"] = option_set_group_assignment_dict[option_set] if option_set_group_assignment_dict is not None else None
            if ("purge_dups" in config["qc_settings"]) and ("qc_datatypes" in config["qc_settings"]["purge_dups"]) and config["qc_settings"]["purge_dups"]["qc_datatypes"]:
                stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"] = config["qc_settings"]["purge_dups"]["qc_datatypes"]
            else:

                if not stage_dict["contig"]["parameters"][parameters_label]["option_set"]["qc_datatypes"]:
                    stage_dict["contig"]["parameters"][parameters_label]["option_set"]["qc_datatypes"] = stage_dict["contig"]["parameters"][parameters_label]["option_set"]["main_datatypes"]
                if not stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"]:
                    stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"] = stage_dict[prev_stage]["parameters"][prev_parameters]["option_set"]["main_datatypes"]

if "polishing" in config["stage_list"]:
    current_stage = "polishing"
    #stage_dict[current_stage] = {}
    prev_stage = stage_dict[current_stage]["prev_stage"]
    tool_list = config["stage_coretools"][current_stage]["default"]
    stage_dict[current_stage]["parameters"] = {}



    for tool in tool_list:
        option_set_group_dict, option_set_group_assignment_dict = None, None
        for option_set in config["coretool_option_sets"][tool]:
            for prev_parameters in stage_dict[prev_stage]["parameters"]:
                parameters_label = "{0}..{1}_{2}".format(prev_parameters, tool, option_set)
                stage_dict[current_stage]["parameters"][parameters_label] = {}
                stage_dict[current_stage]["parameters"][parameters_label]["stage_seq_type"] = "scaffold"
                stage_dict[current_stage]["parameters"][parameters_label]["included"] = True
                stage_dict[current_stage]["parameters"][parameters_label]["prev_stage"] = prev_stage
                stage_dict[current_stage]["parameters"][parameters_label]["prev_parameters"] = prev_parameters
                stage_dict[current_stage]["parameters"][parameters_label]["tool"] = tool
                stage_dict[current_stage]["parameters"][parameters_label]["option_set"] = parameters["tool_options"][tool][option_set]
                stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"] = stage_dict[stage_dict[current_stage]["prev_stage"]]["parameters"][prev_parameters]["haplotype_list"]

                if (len(stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]) == 1) and (stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["use_phased_reads"]):
                    #stage_dict["hic_scaffolding"]["parameters"][parameters_label]["included"] = False
                    stage_dict[current_stage]["parameters"].pop(parameters_label)
                    print(f"WARNING!!! Impossible to phase reads for {parameters_label} as input draft assembly is haploid")
                if not stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["qc_datatypes"]:
                    stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["qc_datatypes"] = stage_dict[prev_stage]["parameters"][prev_parameters]["option_set"]["qc_datatypes"]

                if ("purge_dups" in config["qc_settings"]) and ("qc_datatypes" in config["qc_settings"]["purge_dups"]) and config["qc_settings"]["purge_dups"]["qc_datatypes"]:
                    stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"] = config["qc_settings"]["purge_dups"]["qc_datatypes"]
                else:
                    if ("purge_dups_qc_datatypes" not in stage_dict[current_stage]["parameters"][parameters_label]["option_set"]) or (not stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"]):
                        stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"] = stage_dict[prev_stage]["parameters"][prev_parameters]["option_set"]["purge_dups_qc_datatypes"]

    parameters_list = list(stage_dict[current_stage]["parameters"].keys())

    results_list += [
                    *[expand(out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["polishing", ],
                           haplotype=stage_dict["polishing"]["parameters"][parameters_label]["haplotype_list"],
                           parameters=[parameters_label]) for parameters_label in stage_dict["polishing"]["parameters"]],
                    expand(out_dir_path / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["polishing"],),
                    ]

#for option_supergroup in ["options_affecting_error_correction"]:
            #    stage_dict["contig"]["parameters"][parameters_label][option_supergroup] = option_cluster_reverse_dict[assembler][option_supergroup][option_set]
#print (stage_dict)
if "contig" in config["stage_list"]:
    current_stage = "contig"
    parameters_list = list(stage_dict["contig"]["parameters"].keys())
    print(stage_dict["contig"]["parameters"])
    """
    for parameter_label in stage_dict["contig"]["parameters"]:
        print(parameters_label)
        read_filelist = []
        #print(parameters["tool_options"]["hifiasm"][parameters_label])
        for datatype in stage_dict["contig"]["parameters"][parameters_label]["option_set"]["main_datatypes"]:
            if datatype not in input_filedict:
                continue
            read_filelist += expand(
                output_dict["data"] / ("fastq/{datatype}/filtered/{fileprefix}%s" % config["fastq_extension"]),
                fileprefix=input_file_prefix_dict[datatype],
                datatype=[datatype, ],
                allow_missing=True)
        print(read_filelist)
    """
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
    #print([parameters["tool_options"]["assembly_qc"]["gap"]])
    results_list += [
                     [[[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/trackplots/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.{scaffold_length}.win{window}.step{step}.{threshold_type}.png",
                               scaffold_length=config["qc_settings"]["assembly_scaffold_sets"],
                               threshold_type=config["qc_settings"]["threshold_types"],
                               genome_prefix=[config["genome_prefix"], ],
                               assembly_stage=[current_stage, ],
                               track_type=[track_type],
                               window=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["window"]],
                               step=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["step"]],
                               haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                               parameters=[parameters_label])

                        for window_settings in config["qc_settings"]["windows_sets"]]
                        for parameters_label in stage_dict[current_stage]["parameters"]]
                        for track_type in ["gc"]],  #"windowmasker", "trf"
                     [expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.bedgraph",
                            genome_prefix=[config["genome_prefix"], ],
                            assembly_stage=[current_stage, ],
                            parameters=[parameters_label],
                            haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                           ) for parameters_label in stage_dict[current_stage]["parameters"]
                     ]
                    ]
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
                            *[expand(out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs_adaptor/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.unfiltered.{database}.report",
                                   genome_prefix=[config["genome_prefix"], ],
                                   assembly_stage=["contig"],
                                   haplotype=stage_dict["contig"]["parameters"][parameters_label]["haplotype_list"],
                                   parameters=[parameters_label],
                                   database=config["database_set"]["fcs_adaptor"]) for parameters_label in parameters_list],
                            ]
        if config["database_set"]["fcs"] and (not config["skip_fcs"]):
            results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/contamination_scan/{haplotype}/fcs/{database}/{genome_prefix}.{assembly_stage}.{haplotype}.unfiltered.{database}.taxonomy",
                                    genome_prefix=[config["genome_prefix"], ],
                                    assembly_stage=["contig"],
                                    haplotype=stage_dict["contig"]["parameters"][parameters_label]["haplotype_list"] + (["alt" if stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] > 1 else "alt0"] if stage_dict["contig"]["parameters"][parameters_label]["assembler"] == "hifiasm" else []), # TODO: modify "alt" when assemblers other than hifiasm will be added
                                    parameters=[parameters_label],
                                    database=config["database_set"]["fcs"]) for parameters_label in parameters_list]
                            ]
    if (not config["skip_combined_hic"]) and (not config["skip_combined_contig_hic"]) and ("hic" in data_types):
        results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.NA.{haplotype}.rmdup.pre.mapq{mapq}.hic",
                                      assembly_stage=["contig"],
                                      parameters=[parameters_label],
                                      genome_prefix=[config["genome_prefix"], ],
                                      haplotype=["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                                mapq=parameters["tool_options"]["yahs_juicer_pre"]["mapq"]) if stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] > 1 else []
                               for parameters_label in parameters_list] if not config["skip_hic_file"] else []]
    if current_stage in config["extended_qc_stages"]:
        results_list += [[[
                              expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.NA.{haplotype}.{subset}.rmdup.mapq{mapq}.{res}.tracks.pretext",
                                  res=["high_res"], #"default",
                                  haplotype=["reordered" if ("bird_genome" in config) and config[
                                      "bird_genome"] else "combined"],
                                  subset=["all"] , #+ (["microchr"] if ("bird_genome" in config) and config["bird_genome"] else [])
                                  genome_prefix=[config["genome_prefix"], ],
                                  assembly_stage=[current_stage],
                                  parameters=stage_dict[current_stage]["parameters"],
                                  resolution=parameters["tool_options"]["pretextsnapshot"]["resolution"],
                                  mapq=parameters["tool_options"]["pretextmap"]["mapq"],
                                  ext=parameters["tool_options"]["pretextsnapshot"]["format"],
                                  #window=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["window"],
                                  #step = parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["step"],
                                     ) for window_step_set in config["qc_settings"]["windows_sets"]] for
                          parameters_label in
                          stage_dict[current_stage]["parameters"]] if coverage_track_data_type_set else [],
                         ]

if "purge_dups" in config["stage_list"]:
    current_stage = "purge_dups"
    prev_stage = stage_dict[current_stage]["prev_stage"]
    purge_dupser_list = config["stage_coretools"][current_stage]["default"]
    stage_dict[current_stage]["parameters"] = {}
    for purge_dupser in purge_dupser_list:
        for option_set in config["coretool_option_sets"][purge_dupser]:
            for prev_parameters in stage_dict[prev_stage]["parameters"]:
                parameters_label = "{0}..{1}_{2}".format(prev_parameters, purge_dupser, option_set)
                stage_dict[current_stage]["parameters"][parameters_label] = {}
                stage_dict[current_stage]["parameters"][parameters_label]["included"] = True
                stage_dict[current_stage]["parameters"][parameters_label]["prev_stage"] = prev_stage
                stage_dict[current_stage]["parameters"][parameters_label]["prev_parameters"] = prev_parameters
                stage_dict[current_stage]["parameters"][parameters_label]["purge_dupser"] = purge_dupser
                stage_dict[current_stage]["parameters"][parameters_label]["stage_seq_type"] = "contig"
                stage_dict[current_stage]["parameters"][parameters_label]["option_set"] = parameters["tool_options"][purge_dupser][option_set]
                stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"] = stage_dict[stage_dict[current_stage]["prev_stage"]]["parameters"][prev_parameters]["haplotype_list"]
                if not stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["main_datatypes"]:
                    stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["main_datatypes"] = stage_dict[stage_dict[current_stage]["prev_stage"]]["parameters"][prev_parameters]["option_set"]["main_datatypes"]
                if not stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["qc_datatypes"]:
                    stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["qc_datatypes"] = stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["main_datatypes"]

                if ("purge_dups" in config["qc_settings"]) and ("qc_datatypes" in config["qc_settings"]["purge_dups"]) and config["qc_settings"]["purge_dups"]["qc_datatypes"]:
                    stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"] = config["qc_settings"]["purge_dups"]["qc_datatypes"]
                else:
                    if ("purge_dups_qc_datatypes" not in stage_dict[current_stage]["parameters"][parameters_label]["option_set"]) or (not stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"]):
                        stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"] = stage_dict[prev_stage]["parameters"][prev_parameters]["option_set"]["purge_dups_qc_datatypes"]

    parameters_list = list(stage_dict[current_stage]["parameters"].keys())
    results_list += [
                     *[expand(out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.purge_dups.{haplotype}.fasta",
                              genome_prefix=[config["genome_prefix"], ],
                              assembly_stage=[current_stage],
                              haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                              parameters=[parameters_label]) for parameters_label in parameters_list],
                    *[expand(out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
                             genome_prefix=[config["genome_prefix"], ],
                             assembly_stage=[current_stage],
                             haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                             parameters=[parameters_label]) for parameters_label in parameters_list],
                    *[expand(out_dir_path /  "{assembly_stage}/{parameters}/assembly_qc/purge_dups/{haplotype}/PB.stat",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=[current_stage],
                           haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                           parameters=[parameters_label]) for parameters_label in parameters_list],
                    expand(out_dir_path / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=[current_stage],),
                    expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/purge_dups/after.comparison.coverage.png",
                        assembly_stage=[current_stage],
                        parameters=parameters_list
                           ),

                    ]
    results_list += [
                     [[[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/trackplots/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.{scaffold_length}.win{window}.step{step}.{threshold_type}.png",
                               scaffold_length=config["qc_settings"]["assembly_scaffold_sets"],
                               threshold_type=config["qc_settings"]["threshold_types"],
                               genome_prefix=[config["genome_prefix"], ],
                               assembly_stage=[current_stage, ],
                               track_type=[track_type],
                               window=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["window"]],
                               step=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["step"]],
                               haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                               parameters=[parameters_label])

                        for window_settings in config["qc_settings"]["windows_sets"]]
                        for parameters_label in stage_dict[current_stage]["parameters"]]
                        for track_type in ["gc"]],  #"windowmasker", "trf"
                     [expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.bedgraph",
                            genome_prefix=[config["genome_prefix"], ],
                            assembly_stage=[current_stage, ],
                            parameters=[parameters_label],
                            haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                           ) for parameters_label in stage_dict[current_stage]["parameters"]
                     ]
                    ]

    if current_stage in config["extended_qc_stages"]:
        results_list += [[[[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/trackplots/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.{scaffold_length}.win{window}.step{step}.{threshold_type}.png",
                               scaffold_length=config["qc_settings"]["assembly_scaffold_sets"],
                               threshold_type=config["qc_settings"]["threshold_types"],
                               genome_prefix=[config["genome_prefix"], ],
                               assembly_stage=[current_stage, ],
                               track_type=[track_type],
                               window=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["window"]],
                               step=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["step"]],
                               haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                               parameters=[parameters_label])

                        for window_settings in config["qc_settings"]["windows_sets"]]
                        for parameters_label in stage_dict[current_stage]["parameters"]]
                        for track_type in ["windowmasker"] + (["trf"] if not config["skip_trf"] else [])],]
        if not config["skip_wga"]:
            results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/wga.{query_prefix}.{query_length}.to.{target_prefix}.{target_length}.YASS.R11.soft.min_len{min_target_len}.png",
                                     query_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     target_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     genome_prefix=[config["genome_prefix"], ],
                                     assembly_stage=[current_stage, ],
                                     parameters=[parameters_label],
                                     min_target_len=parameters["tool_options"]["wga"]["min_target_len"],
                                     query_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                     target_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                   ) for parameters_label in stage_dict[current_stage]["parameters"]],
                             [expand(out_dir_path / "{assembly_stage}/{parameters}/wga.{query_prefix}.{query_length}.to.{target_prefix}.{target_length}.YASS.R11.soft.min_len{min_target_len}.png",
                                     query_length=config["qc_settings"]["reference_scaffold_sets"],
                                     target_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     genome_prefix=[config["genome_prefix"], ],
                                     assembly_stage=[current_stage, ],
                                     parameters=[parameters_label],
                                     min_target_len=parameters["tool_options"]["wga"]["min_target_len"],
                                     query_prefix=list(input_reference_filedict.keys()),
                                     target_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                   ) for parameters_label in stage_dict[current_stage]["parameters"]],
                             ]

        if coverage_track_data_type_set:
            results_list += [[[
                                  expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/trackplots/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.coverage_{settings}.{scaffold_length}.win{window}.step{step}.png",
                                      settings=parameters["tool_options"]["mosdepth"]["options"],
                                      scaffold_length=config["qc_settings"]["assembly_scaffold_sets"],
                                      window=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["window"],
                                      step=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["step"],
                                      genome_prefix=[config["genome_prefix"], ],
                                      assembly_stage=[current_stage, ],
                                      datatype=coverage_track_data_type_set,
                                      haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                                      parameters=[parameters_label]) for window_step_set in
                                            config["qc_settings"]["windows_sets"]] for parameters_label in
                                                      stage_dict[current_stage]["parameters"]] if coverage_track_data_type_set else [],
                             ]

    for parameters_label in parameters_list:
        if "skipped" not in parameters_label:
            results_list += [[expand(out_dir_path / "purge_dups/{parameters}/{purge_stage}/{haplotype}/{genome_prefix}.dups.{artefact}.fasta",
                                    purge_stage=["first_stage", ] if haplotype == "hap0" else ["first_stage", "second_stage"],
                                    genome_prefix=[config["genome_prefix"], ],
                                    artefact=["junk", "repeat", "haplotig", "ovlp", "highcov"],
                                    haplotype=[haplotype],
                                    parameters=[parameters_label]) for haplotype in stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"]],
                             expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/purge_dups/before.comparison.coverage.png",
                                 assembly_stage=["purge_dups"],
                                 parameters=[parameters_label]
                                    ),
                             [expand(out_dir_path /  "purge_dups/{parameters}/assembly_qc/purge_dups/{haplotype}/{haplotype}.before-after.comparison.coverage.png",
                                    parameters=[parameters_label],
                                    haplotype=[haplotype],
                                    ) for haplotype in stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"]]
                             ]

    if not config["skip_busco"]:
        results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.busco5.{busco_lineage}.tar.gz",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=[current_stage, ],
                                haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=[current_stage],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=[current_stage],
                                haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=[current_stage],
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
    current_stage = "hic_scaffolding"
    prev_stage = stage_dict[current_stage]["prev_stage"]
    hic_scaffolder_list = config["stage_coretools"][current_stage]["default"]
    stage_dict[current_stage]["parameters"] = {}

    for hic_scaffolder in hic_scaffolder_list:
        for option_set in config["coretool_option_sets"][hic_scaffolder]:
            for prev_parameters in stage_dict[prev_stage]["parameters"]:
                parameters_label = "{0}..{1}_{2}".format(prev_parameters, hic_scaffolder, option_set)
                stage_dict[current_stage]["parameters"][parameters_label] = {}
                stage_dict[current_stage]["parameters"][parameters_label]["stage_seq_type"] = "scaffold"
                stage_dict[current_stage]["parameters"][parameters_label]["included"] = True
                stage_dict[current_stage]["parameters"][parameters_label]["prev_stage"] = prev_stage
                stage_dict[current_stage]["parameters"][parameters_label]["prev_parameters"] = prev_parameters
                stage_dict[current_stage]["parameters"][parameters_label]["hic_scaffolder"] = hic_scaffolder
                stage_dict[current_stage]["parameters"][parameters_label]["option_set"] = parameters["tool_options"][hic_scaffolder][option_set]
                stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"] = stage_dict[stage_dict[current_stage]["prev_stage"]]["parameters"][prev_parameters]["haplotype_list"]

                if (len(stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]) == 1) and (stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["use_phased_reads"]):
                    #stage_dict["hic_scaffolding"]["parameters"][parameters_label]["included"] = False
                    stage_dict[current_stage]["parameters"].pop(parameters_label)
                    print(f"WARNING!!! Impossible to phase reads for {parameters_label} as input draft assembly is haploid")
                if not stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["qc_datatypes"]:
                    stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["qc_datatypes"] = stage_dict[prev_stage]["parameters"][prev_parameters]["option_set"]["qc_datatypes"]

                if ("purge_dups" in config["qc_settings"]) and ("qc_datatypes" in config["qc_settings"]["purge_dups"]) and config["qc_settings"]["purge_dups"]["qc_datatypes"]:
                    stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"] = config["qc_settings"]["purge_dups"]["qc_datatypes"]
                else:
                    if ("purge_dups_qc_datatypes" not in stage_dict[current_stage]["parameters"][parameters_label]["option_set"]) or (not stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"]):
                        stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"] = stage_dict[prev_stage]["parameters"][prev_parameters]["option_set"]["purge_dups_qc_datatypes"]


    #for parameter_label in stage_dict["hic_scaffolding"]["parameters"].keys(): # remove ignore
    #    if not stage_dict["hic_scaffolding"]["parameters"][parameter_label]["included"]:
    #        stage_dict["hic_scaffolding"]["parameters"].pop(parameter_label)

    parameters_list = list(stage_dict[current_stage]["parameters"].keys())

    #haplotype=stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"],
    #                            parameters=[parameters_label]) for parameters_label in parameters_list]
    # stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"]
    #prev_parameters_label = stage_dict["hic_scaffolding"]["parameters"][parameters_label]["prev_parameters"]

    if config["other_tool_option_sets"]["mapping_pipeline"] != "juicer":
        results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam.general_stats",
                            genome_prefix=[config["genome_prefix"], ],
                            assembly_stage=[prev_stage,],
                            haplotype=stage_dict[prev_stage]["parameters"][stage_dict["hic_scaffolding"]["parameters"][current_parameter_label]["prev_parameters"]]["haplotype_list"],
                            phasing_kmer_length=[stage_dict["hic_scaffolding"]["parameters"][current_parameter_label]["option_set"]["phasing_kmer_length"]], #[stage_dict["hic_scaffolding"]["parameters"][parameters_label]["option_set"]["phasing_kmer_length"] for parameter_label in stage_dict["hic_scaffolding"]["parameters"]],
                            parameters=[stage_dict["hic_scaffolding"]["parameters"][current_parameter_label]["prev_parameters"]],) if "threeddna" not in current_parameter_label else [] for current_parameter_label in stage_dict["hic_scaffolding"]["parameters"]],
                         ]

    results_list += [
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

    results_list += [
                     [[[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/trackplots/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.{scaffold_length}.win{window}.step{step}.{threshold_type}.png",
                               scaffold_length=config["qc_settings"]["assembly_scaffold_sets"],
                               threshold_type=config["qc_settings"]["threshold_types"],
                               genome_prefix=[config["genome_prefix"], ],
                               assembly_stage=[current_stage, ],
                               track_type=[track_type],
                               window=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["window"]],
                               step=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["step"]],
                               haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"] + ["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                               parameters=[parameters_label])

                        for window_settings in config["qc_settings"]["windows_sets"]]
                        for parameters_label in stage_dict[current_stage]["parameters"]]
                        for track_type in ["gc"]],  #"windowmasker", "trf"
                     [expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.bedgraph",
                            genome_prefix=[config["genome_prefix"], ],
                            assembly_stage=[current_stage, ],
                            parameters=[parameters_label],
                            haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"] + ["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                           ) for parameters_label in stage_dict[current_stage]["parameters"]
                     ]
                    ]
    if not config["skip_higlass"]:
        results_list += [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.NA.{haplotype}.nodup.pairs",
                                haplotype=["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=[current_stage, ],
                                parameters=stage_dict[current_stage]["parameters"],
                                ),
                         expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.NA.{haplotype}.nodup.higlass.mcool",
                                haplotype=["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=[current_stage, ],
                                parameters=stage_dict[current_stage]["parameters"],
                                ),]
    if not config["skip_combined_hic"]:
        results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.NA.{haplotype}.rmdup.pre.mapq{mapq}.hic",
                                      assembly_stage=["hic_scaffolding"],
                                      parameters=[parameters_label],
                                      genome_prefix=[config["genome_prefix"], ],
                                      haplotype=["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                                mapq=parameters["tool_options"]["yahs_juicer_pre"]["mapq"])
                               for parameters_label in parameters_list] if not config["skip_hic_file"] else []]

    if not (config["skip_prescaf_pretext"] or config["skip_both_pretext"]):
        results_list += [*[expand(#out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.mapq{mapq}.{resolution}.{ext}",
                            out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.all.mapq{mapq}.default.{resolution}.{ext}",
                                      genome_prefix=[config["genome_prefix"], ],
                                      assembly_stage=[prev_stage,],
                                      haplotype=stage_dict[prev_stage]["parameters"][stage_dict["hic_scaffolding"]["parameters"][current_parameter_label]["prev_parameters"]]["haplotype_list"],
                                      phasing_kmer_length=[stage_dict["hic_scaffolding"]["parameters"][current_parameter_label]["option_set"]["phasing_kmer_length"]], #[stage_dict["hic_scaffolding"]["parameters"][parameters_label]["option_set"]["phasing_kmer_length"] for parameter_label in stage_dict["hic_scaffolding"]["parameters"]],
                                      parameters=[stage_dict["hic_scaffolding"]["parameters"][current_parameter_label]["prev_parameters"]],
                                      mapq=parameters["tool_options"]["pretextmap"]["mapq"],
                                      resolution=parameters["tool_options"]["pretextsnapshot"]["resolution"],
                                      ext=parameters["tool_options"]["pretextsnapshot"]["format"])  for current_parameter_label in stage_dict["hic_scaffolding"]["parameters"]],
                            ]

    results_list += [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.NA.{haplotype}.{subset}.rmdup.mapq{mapq}.{res}.pretext",
                                  res=["high_res"], #"default",
                                  haplotype=["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                                  subset=["all"] , #+ ( ["microchr"] if ("bird_genome" in config) and config["bird_genome"] else [] )
                                  genome_prefix=[config["genome_prefix"], ],
                                  assembly_stage=[current_stage],
                                  parameters=stage_dict[current_stage]["parameters"],
                                  resolution=parameters["tool_options"]["pretextsnapshot"]["resolution"],
                                  mapq=parameters["tool_options"]["pretextmap"]["mapq"],
                                  ext=parameters["tool_options"]["pretextsnapshot"]["format"]),
                    expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.NA.{haplotype}.{subset}.mapq{mapq}.default.{resolution}.{ext}",
                                  haplotype=["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                                  subset=["all"] , # + ( ["microchr"] if ("bird_genome" in config) and config["bird_genome"] else [] )
                                  genome_prefix=[config["genome_prefix"], ],
                                  assembly_stage=[current_stage],
                                  parameters=stage_dict[current_stage]["parameters"],
                                  resolution=parameters["tool_options"]["pretextsnapshot"]["resolution"],
                                  mapq=parameters["tool_options"]["pretextmap"]["mapq"],
                                  ext=parameters["tool_options"]["pretextsnapshot"]["format"]),

                        ]

    results_list += [[[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/{genome_prefix}.{assembly_stage}.NA.{haplotype}.{subset}.rmdup.mapq{mapq}.{res}.tracks.pretext",
                              res=["high_res"], #"default",
                              haplotype=["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                              subset=["all"] , # + (["microchr"] if ("bird_genome" in config) and config["bird_genome"] else [])
                              genome_prefix=[config["genome_prefix"], ],
                              assembly_stage=[current_stage],
                              parameters=stage_dict[current_stage]["parameters"],
                              resolution=parameters["tool_options"]["pretextsnapshot"]["resolution"],
                              mapq=parameters["tool_options"]["pretextmap"]["mapq"],
                              ext=parameters["tool_options"]["pretextsnapshot"]["format"],
                              #window=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["window"],
                              #step = parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["step"],
                             ) for window_step_set in config["qc_settings"]["windows_sets"]] for parameters_label in
                                                      stage_dict[current_stage]["parameters"]] if coverage_track_data_type_set else [],
                     ]

    if current_stage in config["extended_qc_stages"]:
        results_list += [[[[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/trackplots/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.{scaffold_length}.win{window}.step{step}.{threshold_type}.png",
                               scaffold_length=config["qc_settings"]["assembly_scaffold_sets"],
                               threshold_type=config["qc_settings"]["threshold_types"],
                               genome_prefix=[config["genome_prefix"], ],
                               assembly_stage=[current_stage, ],
                               track_type=[track_type],
                               window=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["window"]],
                               step=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["step"]],
                               haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"] + ["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                               parameters=[parameters_label])

                        for window_settings in config["qc_settings"]["windows_sets"]]
                        for parameters_label in stage_dict[current_stage]["parameters"]]
                        for track_type in ["windowmasker"] + (["trf"] if not config["skip_trf"] else [])],]
        if not config["skip_purge_dups_qc"]:
            print(stage_dict[current_stage]["parameters"][parameters_label]["option_set"])
            print(data_types)
            results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/purge_dups/{haplotype}/{datatype}/{genome_prefix}.{assembly_stage}.{haplotype}.dups.extended.bed",
                                    assembly_stage=[current_stage,],
                                    parameters=[parameters_label],
                                    haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                                    datatype=set(stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["purge_dups_qc_datatypes"]) & set(data_types) ,
                                    genome_prefix=[config["genome_prefix"], ],
                                    ) for parameters_label in stage_dict[current_stage]["parameters"]]]


        if not config["skip_wga"]:
            results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/wga.{query_prefix}.{query_length}.to.{target_prefix}.{target_length}.YASS.R11.soft.min_len{min_target_len}.png",
                                     query_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     target_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     genome_prefix=[config["genome_prefix"], ],
                                     assembly_stage=[current_stage, ],
                                     parameters=[parameters_label],
                                     min_target_len=parameters["tool_options"]["wga"]["min_target_len"],
                                     query_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                     target_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                   ) for parameters_label in stage_dict[current_stage]["parameters"]],
                             [expand(out_dir_path / "{assembly_stage}/{parameters}/wga.{query_prefix}.{query_length}.to.{target_prefix}.{target_length}.YASS.R11.soft.min_len{min_target_len}.png",
                                     query_length=config["qc_settings"]["reference_scaffold_sets"],
                                     target_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     genome_prefix=[config["genome_prefix"], ],
                                     assembly_stage=[current_stage, ],
                                     parameters=[parameters_label],
                                     min_target_len=parameters["tool_options"]["wga"]["min_target_len"],
                                     query_prefix=list(input_reference_filedict.keys()),
                                     target_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                   ) for parameters_label in stage_dict[current_stage]["parameters"]],
                             ]

        if coverage_track_data_type_set:
            results_list += [[[
                                  expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/trackplots/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.coverage_{settings}.{scaffold_length}.win{window}.step{step}.png",
                                      scaffold_length=config["qc_settings"]["assembly_scaffold_sets"],
                                      window=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["window"],
                                      step=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["step"],
                                      genome_prefix=[config["genome_prefix"], ],
                                      assembly_stage=[current_stage, ],
                                      datatype=coverage_track_data_type_set,
                                      settings=parameters["tool_options"]["mosdepth"]["options"],
                                      haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"] + ["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                                      parameters=[parameters_label]) for window_step_set in
                                            config["qc_settings"]["windows_sets"]] for parameters_label in
                                                      stage_dict[current_stage]["parameters"]] if coverage_track_data_type_set else [],
                             ]

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

    if ("bird_genome" in config) and config["bird_genome"]:
        #print(current_stage)
        results_list += [[expand(out_dir_path / "hic_scaffolding/{parameters}/{genome_prefix}.{assembly_stage}.candidates.microchromosomes.filtered.tsv",
                                 assembly_stage=["hic_scaffolding"],
                                 parameters=[parameter_label],
                                 #haplotype=["combined"],
                                 genome_prefix=[config["genome_prefix"],],
                                 #max_length=[parameters["tool_options"]["microsome_detection"]["max_length"]],
                                 ) for parameter_label in stage_dict["hic_scaffolding"]["parameters"]]
                         ]

    if candidate_agp_filename:
        results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/NA/per_chr/{genome_prefix}.{assembly_stage}.{haplotype}.NA.{candidate_chr_id}.rmdup.precurated.mapq{mapq}.{res}.tracks.pretext",
                                candidate_chr_id=candidate_chr_id_list,
                                assembly_stage=["hic_scaffolding"],
                                parameters=[parameter_label],
                                haplotype=["reordered" if ("bird_genome" in config) and config["bird_genome"] else "combined"],
                                genome_prefix=[config["genome_prefix"], ],
                                res=["high_res"],
                                mapq=parameters["tool_options"]["pretextmap"]["mapq"],
                                #window=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["window"],
                                #step=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["step"],
                                ) for parameter_label in stage_dict[current_stage]["parameters"] for window_step_set in config["qc_settings"]["windows_sets"]]
                         ]



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

if "ref_scaffolding" in config["stage_list"]:
    current_stage = "ref_scaffolding"
    prev_stage = stage_dict[current_stage]["prev_stage"]
    curation_tool_list = config["stage_coretools"][current_stage]["default"]
    stage_dict[current_stage]["parameters"] = {}

    ref_scaffolding_tool_list = config["stage_coretools"][current_stage]["default"]
    stage_dict[current_stage]["parameters"] = {}
    for ref_scaffolding_tool in ref_scaffolding_tool_list:
        for option_set in config["coretool_option_sets"][ref_scaffolding_tool]:
            #print(prev_stage)
            for prev_parameters in stage_dict[prev_stage]["parameters"]:
                for reference in list(input_reference_filedict.keys()):
                    parameters_label = "{0}..{1}_{2}@{3}".format(prev_parameters, ref_scaffolding_tool, option_set, reference)
                    stage_dict[current_stage]["parameters"][parameters_label] = {}
                    stage_dict[current_stage]["parameters"][parameters_label]["stage_seq_type"] = "scaffold"
                    stage_dict[current_stage]["parameters"][parameters_label]["included"] = True
                    stage_dict[current_stage]["parameters"][parameters_label]["ref_scaffolder"] = ref_scaffolding_tool
                    stage_dict[current_stage]["parameters"][parameters_label]["prev_stage"] = prev_stage
                    stage_dict[current_stage]["parameters"][parameters_label]["prev_parameters"] = prev_parameters
                    stage_dict[current_stage]["parameters"][parameters_label]["option_set"] = parameters["tool_options"][ref_scaffolding_tool][option_set] if ref_scaffolding_tool in parameters["tool_options"] else None
                    stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"] = stage_dict[stage_dict[current_stage]["prev_stage"]]["parameters"][prev_parameters]["haplotype_list"]
                    if not stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["qc_datatypes"]:
                        stage_dict[current_stage]["parameters"][parameters_label]["option_set"]["qc_datatypes"] = stage_dict[prev_stage]["parameters"][prev_parameters]["option_set"]["qc_datatypes"]

    parameters_list = list(stage_dict[current_stage]["parameters"].keys())
    #print(stage_dict["ref_scaffolding"])
    results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",#"{assembly_stage}/{parameters}/{genome_prefix}.ref_scaffolding.{haplotype}.fasta",
                             assembly_stage=[current_stage],
                             parameters=[parameters_label],
                             genome_prefix=[config["genome_prefix"], ],
                             haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                             ) for parameters_label in parameters_list],
                     [expand(out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=[current_stage],
                                haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                     expand(out_dir_path / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats",
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=[current_stage],),]

    #results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{genome_prefix}.ref_scaffolding.{haplotype}.fasta",
    #                        genome_prefix=[config["genome_prefix"], ],
    #                        assembly_stage=[current_stage],
    #                        haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
    #                        parameters=[parameters_label]) for parameters_label in parameters_list
    #
    #                  ]]

    #print(parameters["tool_options"]["assembly_qc"])

    #out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.win{window}.step{step}.track.bedgraph",
    #out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/track_stats/{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.win{window}.step{step}.track.stat",
    #out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/track_stats/{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.win{window}.step{step}.track.thresholds",

    results_list += [
                     [[[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/trackplots/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.{scaffold_length}.win{window}.step{step}.{threshold_type}.png",
                               scaffold_length=config["qc_settings"]["assembly_scaffold_sets"],
                               threshold_type=config["qc_settings"]["threshold_types"],
                               genome_prefix=[config["genome_prefix"], ],
                               assembly_stage=[current_stage, ],
                               track_type=[track_type],
                               window=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["window"]],
                               step=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["step"]],
                               haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                               parameters=[parameters_label])

                        for window_settings in config["qc_settings"]["windows_sets"]]
                        for parameters_label in stage_dict[current_stage]["parameters"]]
                        for track_type in ["gc"]],  #"windowmasker", "trf"
                     [expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.canonical_telomere.win1000.step200.track.bedgraph",
                            genome_prefix=[config["genome_prefix"], ],
                            assembly_stage=[current_stage, ],
                            parameters=[parameters_label],
                            haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                           ) for parameters_label in stage_dict[current_stage]["parameters"]
                     ]
                    ]
    if current_stage in config["extended_qc_stages"]:
        results_list += [[[[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/trackplots/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.{scaffold_length}.win{window}.step{step}.{threshold_type}.png",
                               scaffold_length=config["qc_settings"]["assembly_scaffold_sets"],
                               threshold_type=config["qc_settings"]["threshold_types"],
                               genome_prefix=[config["genome_prefix"], ],
                               assembly_stage=[current_stage, ],
                               track_type=[track_type],
                               window=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["window"]],
                               step=[parameters["tool_options"]["assembly_qc"][track_type]["options"][window_settings]["step"]],
                               haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                               parameters=[parameters_label])

                        for window_settings in config["qc_settings"]["windows_sets"]]
                        for parameters_label in stage_dict[current_stage]["parameters"]]
                        for track_type in ["windowmasker"] + (["trf"] if not config["skip_trf"] else [])],]
        if not config["skip_wga"]:
            results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/wga.{query_prefix}.{query_length}.to.{target_prefix}.{target_length}.YASS.R11.soft.min_len{min_target_len}.png",
                                     query_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     target_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     genome_prefix=[config["genome_prefix"], ],
                                     assembly_stage=[current_stage, ],
                                     parameters=[parameters_label],
                                     min_target_len=parameters["tool_options"]["wga"]["min_target_len"],
                                     query_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                     target_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                   ) for parameters_label in stage_dict[current_stage]["parameters"]],
                             [expand(out_dir_path / "{assembly_stage}/{parameters}/wga.{query_prefix}.{query_length}.to.{target_prefix}.{target_length}.YASS.R11.soft.min_len{min_target_len}.png",
                                     query_length=config["qc_settings"]["reference_scaffold_sets"],
                                     target_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     genome_prefix=[config["genome_prefix"], ],
                                     assembly_stage=[current_stage, ],
                                     parameters=[parameters_label],
                                     min_target_len=parameters["tool_options"]["wga"]["min_target_len"],
                                     query_prefix=list(input_reference_filedict.keys()),
                                     target_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                   ) for parameters_label in stage_dict[current_stage]["parameters"]],
                             ]

        if coverage_track_data_type_set:
            results_list += [[[
                                  expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/trackplots/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.coverage.{scaffold_length}.win{window}.step{step}.png",
                                      scaffold_length=config["qc_settings"]["assembly_scaffold_sets"],
                                      window=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["window"],
                                      step=parameters["tool_options"]["assembly_qc"]["coverage"]["options"][window_step_set]["step"],
                                      genome_prefix=[config["genome_prefix"], ],
                                      assembly_stage=[current_stage, ],
                                      datatype=coverage_track_data_type_set,
                                      haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"],
                                      parameters=[parameters_label]) for window_step_set in config["qc_settings"]["windows_sets"]]
                                                 for parameters_label in stage_dict[current_stage]["parameters"]] if coverage_track_data_type_set else [],
                             ]


    if not config["skip_busco"]:
        results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.busco5.{busco_lineage}.tar.gz",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["ref_scaffolding", ],
                                haplotype=stage_dict["ref_scaffolding"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                        *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["ref_scaffolding"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                        *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["ref_scaffolding"],
                                haplotype=stage_dict["ref_scaffolding"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                        expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.busco.merged.tsv",
                            busco_lineage=config["busco_lineage_list"],
                            genome_prefix=[config["genome_prefix"], ],
                            assembly_stage=["ref_scaffolding"],
                            parameters=parameters_list
                               ),
                        ]

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
        if not config["skip_gathering"]:
            results_list += [[expand(out_dir_path / "curation_files/{parameters}/{haplotype}/scaffolds",
                                             parameters=[parameter_label],
                                             haplotype=stage_dict["curation"]["parameters"][parameter_label]["haplotype_list"],
                                            ) for parameter_label in stage_dict["curation"]["parameters"]],
                                     ]
        if input_reference_filedict:
            if not config["skip_ragtag"]:
                results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/scaffolds/ragtag/{reference}/{genome_prefix}.{haplotype}.to.{reference}.fasta",
                                        genome_prefix=[config["genome_prefix"], ],
                                        assembly_stage=["curation", ],
                                        parameters=[parameters_label],
                                        reference=list(input_reference_filedict.keys()),
                                        haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                        ) for parameters_label in stage_dict["curation"]["parameters"]],
                                 ]
        if not config["skip_wga"]:
            results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/wga.{query_prefix}.{query_length}.to.{target_prefix}.{target_length}.YASS.R11.soft.min_len{min_target_len}.png",
                                     query_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     target_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     genome_prefix=[config["genome_prefix"], ],
                                     assembly_stage=[current_stage, ],
                                     parameters=[parameters_label],
                                     min_target_len=parameters["tool_options"]["wga"]["min_target_len"],
                                     query_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                     target_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                   ) for parameters_label in stage_dict[current_stage]["parameters"]],
                             [expand(out_dir_path / "{assembly_stage}/{parameters}/wga.{query_prefix}.{query_length}.to.{target_prefix}.{target_length}.YASS.R11.soft.min_len{min_target_len}.png",
                                     query_length=config["qc_settings"]["reference_scaffold_sets"],
                                     target_length=config["qc_settings"]["assembly_scaffold_sets"],
                                     genome_prefix=[config["genome_prefix"], ],
                                     assembly_stage=[current_stage, ],
                                     parameters=[parameters_label],
                                     min_target_len=parameters["tool_options"]["wga"]["min_target_len"],
                                     query_prefix=list(input_reference_filedict.keys()),
                                     target_prefix=expand("{genome_prefix}.{assembly_stage}.{haplotype}",
                                                         genome_prefix=[config["genome_prefix"], ],
                                                         assembly_stage=[current_stage, ],
                                                         haplotype=stage_dict[current_stage]["parameters"][parameters_label]["haplotype_list"]),
                                   ) for parameters_label in stage_dict[current_stage]["parameters"]],
                             ]
        results_list += [[[[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/scaffolds/{genome_prefix}.input.{haplotype}.{track_type}.win{window}.step{step}.{threshold_type}.png",
                                threshold_type=config["qc_settings"]["threshold_types"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                track_type=[track_type],
                                window=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["window"]],
                                step=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["step"]],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for window_settings in stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"] ] for parameters_label in stage_dict["curation"]["parameters"]] for track_type in  ["windowmasker", "gc"] + (["trf"] if not config["skip_trf"] else []) ],
                         [[[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{track_type}.win{window}.step{step}.track.stat",
                                seq_type=["scaffolds"],
                                threshold_type=config["qc_settings"]["threshold_types"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                track_type=[track_type],
                                window=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["window"]],
                                step=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["step"]],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for window_settings in stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"] ] for parameters_label in stage_dict["curation"]["parameters"]] for track_type in ["windowmasker", "gc"] + (["trf"] if not config["skip_trf"] else []) ],
                         [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.canonical.txt",
                                seq_type=["scaffolds"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],
                         [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.canonical_telomere.win1000.step200.track.bedgraph",
                                seq_type=["scaffolds"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],
                         [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.canonical_telomere_warning.win1000.step200.track.bedgraph",
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
        if coverage_track_data_type_set:

            results_list += [[[expand(out_dir_path / "curation/{parameters}/{haplotype}/scaffolds/{genome_prefix}.input.{haplotype}.{datatype}.coverage.win{window}.step{step}.png",
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
                         ]
        if "hic_scaffolding" in config["stage_list"]:
            results_list += [[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.canonical_telomere_warning.win1000.step200.scaled.track.bedgraph",
                                seq_type=["scaffolds"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],
                            [expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.canonical_telomere.win1000.step200.scaled.track.bedgraph",
                                seq_type=["scaffolds"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],
                            [[[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{track_type}.win{window}.step{step}.scaled.track.bedgraph",
                                seq_type=["scaffolds"],
                                threshold_type=config["qc_settings"]["threshold_types"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                track_type=[track_type],
                                window=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["window"]],
                                step=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["step"]],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for window_settings in stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"] ] for parameters_label in stage_dict["curation"]["parameters"]] for track_type in ["windowmasker", "gc"] + (["trf"] if not config["skip_trf"] else []) ],
                             ]
            if coverage_track_data_type_set:
                results_list += [[[expand(out_dir_path / "curation/{parameters}/{haplotype}/scaffolds/{genome_prefix}.input.{haplotype}.{datatype}_{cov_type}_coverage.win{window}.step{step}.scaled.track.bedgraph",
                                cov_type=["mean", "median"],
                                window=stage_dict["curation"]["parameters"][parameters_label]["option_set"]["coverage"]["options"][window_step_set]["window"],
                                step=stage_dict["curation"]["parameters"][parameters_label]["option_set"]["coverage"]["options"][window_step_set]["step"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                datatype=coverage_track_data_type_set,
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for window_step_set in stage_dict["curation"]["parameters"][parameters_label]["option_set"]["coverage"]["options"]] for parameters_label in stage_dict["curation"]["parameters"]] if coverage_track_data_type_set else [],
                            ]

    if ("contigs" in config["curation_seq_type"]) and (prev_stage == "hic_scaffolding") :
        results_list += [[[[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{track_type}.win{window}.step{step}.track.{filetype}",
                                filetype=["stat", "assembly.bedgraph"],
                                seq_type=["contigs"],
                                threshold_type=config["qc_settings"]["threshold_types"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                track_type=[track_type],
                                window=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["window"]],
                                step=[stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"][window_settings]["step"]],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for window_settings in stage_dict["curation"]["parameters"][parameters_label]["option_set"][track_type]["options"] ] for parameters_label in stage_dict["curation"]["parameters"]] for track_type in ["windowmasker", "gc", "gap"] + (["trf"] if not config["skip_trf"] else []) ],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.canonical.txt",
                                seq_type=["contigs"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{telo_type}.win1000.step200.track.assembly.bedgraph",
                                telo_type=["canonical_telomere", "non_canonical_telomere"],
                                seq_type=["contigs"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["curation", ],
                                haplotype=stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in stage_dict["curation"]["parameters"]],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{telo_type}.win1000.step200.track.assembly.bedgraph",
                                telo_type=["canonical_telomere_warning", "non_canonical_telomere_warning"],
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
                                phasing_kmer_length=["31"] ,# [stage_dict["curation"]["parameters"][current_parameter_label]["option_set"]["phasing_kmer_length"]],
                                parameters=[stage_dict["curation"]["parameters"][current_parameter_label]["prev_parameters"]],
                                ) for current_parameter_label in stage_dict["curation"]["parameters"]]

    if (prev_stage == "hic_scaffolding") and (not config["skip_gathering"]):
        #for parameter_label in stage_dict["curation"]["parameters"]:
        #    print(stage_dict["curation"]["parameters"][parameters_label])
        #    print(stage_dict["curation"]["parameters"][parameters_label]["haplotype_list"])
        results_list += [[expand(out_dir_path / "curation_files/{parameters}/{haplotype}/{genome_prefix}.hic_scaffolding.{haplotype}.hic",
                                 parameters=[parameter_label],
                                 haplotype=stage_dict["curation"]["parameters"][parameter_label]["haplotype_list"],
                                 genome_prefix=[config["genome_prefix"], ],
                                 ) for parameter_label in stage_dict["curation"]["parameters"]],
                        [expand(out_dir_path / "curation_files/{parameters}/{haplotype}/contigs",
                                parameters=[parameter_label],
                                haplotype=stage_dict["curation"]["parameters"][parameter_label]["haplotype_list"],
                                ) for parameter_label in stage_dict["curation"]["parameters"]],
                         ]

    if ("bird_genome" in config) and config["bird_genome"]:
        #if "contigs" in config["curation_seq_type"]:
        #    results_list += [[expand(out_dir_path / "curation/{parameters}/{haplotype}/contigs/{genome_prefix}.input.{haplotype}.order.tsv",
        #                             parameters=[parameter_label],
        #                             haplotype=stage_dict["curation"]["parameters"][parameter_label]["haplotype_list"],
        #                             genome_prefix=[config["genome_prefix"], ],
        #                            ) for parameter_label in stage_dict["curation"]["parameters"]],
        #                     ]
        if "scaffolds" in config["curation_seq_type"]:
            results_list += [[expand(out_dir_path / "curation/{parameters}/{haplotype}/scaffolds/{genome_prefix}.input.{haplotype}.max{max_length}.candidates.microchromosomes.filtered.tsv",
                                     parameters=[parameter_label],
                                     haplotype=stage_dict["curation"]["parameters"][parameter_label]["haplotype_list"],
                                     genome_prefix=[config["genome_prefix"],],
                                     max_length=parameters["tool_options"]["microsome_detection"]["max_length"],
                                    ) for parameter_label in stage_dict["curation"]["parameters"]],
                             ]

    with open("tmp.results_list", "w") as out_fd:
        for filename in results_list:
            out_fd.write(str(filename) + "\n")
#----


#---- Final rule ----
pd.Series(results_list).to_csv(config["out_dir"] + "/requested_files.tab", sep="\t", header=False, index=False)
rule all:
    input:
        results_list
        #results_dict[config["mode"]]
#----

#---- Include section ----
include: "workflow/rules/General/Log.smk"
#include: "workflow/rules/Install/Pip.smk"
include: "workflow/rules/Preprocessing/Files.smk"
include: "workflow/rules/Preprocessing/Combine.smk"
include: "workflow/rules/QCFiltering/FastQC.smk"
include: "workflow/rules/QCFiltering/MultiQC.smk"
include: "workflow/rules/QCFiltering/Cutadapt.smk"
include: "workflow/rules/QCFiltering/Trimmomatic.smk"

include: "workflow/rules/QCFiltering/TADbit.smk"
include: "workflow/rules/QCFiltering/Nanopore.smk"
include: "workflow/rules/QCFiltering/NanoQC.smk"
include: "workflow/rules/QCFiltering/NanoPlot.smk"

include: "workflow/rules/Kmer/Jellyfish.smk"
include: "workflow/rules/Kmer/Meryl.smk"
include: "workflow/rules/Kmer/Yak.smk"
include: "workflow/rules/Kmer/Smudgeplot.smk"
include: "workflow/rules/Kmer/GCplot.smk"
include: "workflow/rules/Kmer/Genomescope.smk"
include: "workflow/rules/Kmer/Krater.smk"

include: "workflow/rules/QCAssembly/BUSCO5.smk"
include: "workflow/rules/QCAssembly/Merqury.smk"
include: "workflow/rules/QCAssembly/QUAST.smk"
include: "workflow/rules/QCAssembly/General.smk"
include: "workflow/rules/QCAssembly/Purge_dups.smk"
include: "workflow/rules/Contamination/FCS.smk"
include: "workflow/rules/Contamination/Kraken2.smk"

#if "hifi" in data_types:
include: "workflow/rules/Contigs/Hifiasm.smk"

include: "workflow/rules/Contigs/Graph.smk"
include: "workflow/rules/Stats/General.smk"

if "purge_dups" in config["stage_list"]:
    include: "workflow/rules/Purge_dups/Purge_dups.smk"
    include: "workflow/rules/Purge_dups/Purge_dupsQC.smk"

include: "workflow/rules/HiC/ReadPhasing.smk"

include: "workflow/rules/Alignment/Index.smk"
include: "workflow/rules/Alignment/Common.smk"
include: "workflow/rules/Alignment/Stats.smk"
include: "workflow/rules/Alignment/Winnowmap.smk"

if "hic" in data_types:
    if (sum(list(pd.Series(["hic_scaffolding",
                        "gap_closing",
                        "draft_qc", "contig"]).isin(config["stage_list"]))) > 0) :
        if config["other_tool_option_sets"]["mapping_pipeline"] == "arima":
            print("Mapping pipeline: Arima")
            include: "workflow/rules/Alignment/Arima.smk"
        elif config["other_tool_option_sets"]["mapping_pipeline"] == "bwa_only":
            print("Mapping pipeline: BWA only")
            include: "workflow/rules/Alignment/BWAOnly.smk"
        elif config["other_tool_option_sets"]["mapping_pipeline"] == "pairtools":
            print("Mapping pipeline: Pairtools")
            include: "workflow/rules/Alignment/Pairtools.smk"
        include: "workflow/rules/Alignment/PostAlignment.smk"

    if (sum(list(pd.Series(["hic_scaffolding",
                        "gap_closing",
                        "draft_qc", "contig"]).isin(config["stage_list"]))) > 0) :
        include: "workflow/rules/Alignment/Pretext.smk"

    if "hic_scaffolding" in config["stage_list"]:
        include: "workflow/rules/HiC/YAHS.smk"
        include: "workflow/rules/HiC/3DDNA.smk"

include: "workflow/rules/Polishing/NextPolish2.smk"

include: "workflow/rules/QCAssembly/RapidCuration.smk"
include: "workflow/rules/QCAssembly/GapTrack.smk"
include: "workflow/rules/QCAssembly/WindowmaskerTrack.smk"
include: "workflow/rules/QCAssembly/CoverageTrack.smk"
include: "workflow/rules/QCAssembly/TelomereTrack.smk"
include: "workflow/rules/QCAssembly/TelomereTidkTrack.smk"
include: "workflow/rules/QCAssembly/TRFTrack.smk"
include: "workflow/rules/QCAssembly/Masking.smk"
include: "workflow/rules/QCAssembly/GCTrack.smk"
include: "workflow/rules/QCAssembly/WGA.smk"
include: "workflow/rules/QCAssembly/HiC.smk"
include: "workflow/rules/QCAssembly/MicroChromosomes.smk"
include: "workflow/rules/QCAssembly/PretextPerChr.smk"
include: "workflow/rules/QCAssembly/RagTag.smk"
include: "workflow/rules/mtDNA/MitoHiFi.smk"
include: "workflow/rules/mtDNA/Mitoz.smk"
#include: "workflow/rules/QCAssembly/VariantTrack.smk"

if "gap_closing" in config["stage_list"]:
    include: "workflow/rules/Finalization/GapClosing.smk"

if "ref_scaffolding" in config["stage_list"]:
    include: "workflow/rules/RefScaffolding/RagTag.smk"

#include: "workflow/rules/Curation/CurationFiles.smk"
#include: "workflow/rules/Curation/MicroChromosomes.smk"
