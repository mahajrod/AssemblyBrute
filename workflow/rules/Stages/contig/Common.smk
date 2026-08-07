
localrules: extract_lambda_value

def get_coverage_estimator(wildcards):
    coverage_estimator = stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["coverage_estimator"]
    if "genome_assemblers" in config["tool_manually_adjusted_features"]:
        if "coverage_estimator" in config["tool_manually_adjusted_features"]["genome_assemblers"]:
            if config["tool_manually_adjusted_features"]["genome_assemblers"]["coverage_estimator"] is not None:
                coverage_estimator = config["tool_manually_adjusted_features"]["genome_assemblers"]["coverage_estimator"]

    if coverage_estimator not in config["allowed_coverage_estimators"]:
        raise ValueError("ERROR!!! Unrecognized coverage estimator ({0}) for contig assembly hifiasm_{1} !".format(coverage_estimator,
                                                                                                                   wildcards.parameters))
    return coverage_estimator

def get_coverage_estimator_report_filename(wildcards):
    coverage_estimator = get_coverage_estimator(wildcards)
    report_filename = config["out_dir"] / ("kmer/{0}/final/{1}/{2}.{3}.final.{4}.{5}.{1}.parameters".format("_".join(config["final_kmer_datatypes"]),
                                                                                                               coverage_estimator,
                                                                                                               wildcards.genome_prefix,
                                                                                                               "_".join(config["final_kmer_datatypes"]),
                                                                                                               config["final_kmer_length"],
                                                                                                               config["final_kmer_counter"]))

    return report_filename

rule extract_lambda_value:
    input:
        coverage_estimator_report_filename=get_coverage_estimator_report_filename
    output:
        lambda_file=config["out_dir"]/ "contig/{parameters}/{genome_prefix}.lambda",
    log:
        std=config["out_dir"] / "log/extract_lambda_value.{parameters}.{genome_prefix}.log",
        cluster_log=config["out_dir"] / "log/extract_lambda_value.{parameters}.{genome_prefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/extract_lambda_value.{parameters}.{genome_prefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/extract_lambda_value.{parameters}.{genome_prefix}.benchmark.txt"
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("extract_lambda_value"),
        cpus=parameters["threads"]["extract_lambda_value"],
        time=parameters["time"]["extract_lambda_value"],
        mem=parameters["memory_mb"]["extract_lambda_value"],
    threads:
        parameters["threads"]["extract_lambda_value"]
    run:
        with open(log.std, "w") as log_fd, open(output.lambda_file, "w") as out_fd:
            coverage_estimator = get_coverage_estimator(wildcards)
            if "genome_assemblers" in config["tool_manually_adjusted_features"]:
                if "lambda" in config["tool_manually_adjusted_features"]["genome_assemblers"]:
                    if config["tool_manually_adjusted_features"]["genome_assemblers"]["lambda"] is not None:
                        if isinstance(config["tool_manually_adjusted_features"]["genome_assemblers"]["lambda"], Number):
                            lambda_value = config["tool_manually_adjusted_features"]["genome_assemblers"]["lambda"]
                            print("Using a preset lambda value ({0}) for contig assembly hifiasm_{1} ...".format(config["tool_manually_adjusted_features"]["genome_assemblers"]["lambda"],
                                                                                                                 wildcards.parameters))
                            log_fd.write("Using the preset lambda value ({0}) for contig assembly hifiasm_{1} ...\n".format(lambda_value,
                                                                                                                            wildcards.parameters))

                            log_fd.write("Report file:\tIgnored\n")
                            log_fd.write("Lambda:\t%.2f\n" % lambda_value)
                            out_fd.write("%.2f\n" % lambda_value)
                            return lambda_value
                        else:
                            message = "ERROR!!! Preset lambda value is not a number! Check value in contig['tool_manually_adjusted_features']['hifiasm']['lambda'] ..."
                            log_fd.write(message + "\n")
                            raise ValueError(message)

            log_fd.write("Using the {0} as a coverage estimator for contig assembly hifiasm_{1} ...\n".format(coverage_estimator,
                                                                                                            wildcards.parameters))
            report_filename = config["out_dir"] / ("kmer/{0}/final/{1}/{2}.{3}.final.{4}.{5}.{1}.parameters".format("_".join(config["final_kmer_datatypes"]),
                                                                                                                       coverage_estimator,
                                                                                                                       wildcards.genome_prefix,
                                                                                                                       "_".join(config["final_kmer_datatypes"]),
                                                                                                                       config["final_kmer_length"],
                                                                                                                       config["final_kmer_counter"]))
            log_fd.write("Report file:\t%s\n" % report_filename)
            if coverage_estimator == "genomescope":
                with open(report_filename, "r") as in_fd:
                    for line in in_fd:
                        if "Lambda" in line:
                            lambda_value = float(line.strip().split("\t")[1])
                            break
            elif coverage_estimator == "krater":
                with open(report_filename, "r") as in_fd:
                    for line in in_fd:
                        if "Kmer multiplicity at first maximum" in line:
                            lambda_value = float(line.strip().split("\t")[1])
                            break
                if "krater" in config["tool_manually_adjusted_features"]:
                    if not config["tool_manually_adjusted_features"]["krater"]["use_second_peak"]:
                        lambda_value = lambda_value / 2
                else:
                    lambda_value = lambda_value / 2

            else:
                message = "ERROR!!! Unimplemented (but planned) coverage estimator ({0}) for contig assembly hifiasm_{1} !".format(coverage_estimator,
                                                                                                                                   wildcards.parameters)
                log_fd.write(message + "\n")
                raise ValueError(message)
            log_fd.write("Lambda:\t%.2f\n" % lambda_value)
            out_fd.write("%.2f\n" % lambda_value)

#def get_main_read_filelist(wildcards):
#    read_filelist = []
#    for datatype in stage_dict["contig"].parameters[wildcards.parameters]["option_set"]["main_datatypes"]:
#        if datatype not in config["data"]:
#            continue
#        read_filelist += expand(config["out_dir"]/ "data/{datatype}/final/{fileprefix}{extension}",
#                                fileprefix=config["data"][datatype]["conv_file_prefix_list"],
#                                datatype=[datatype,],
#                                extension=[config["data"][datatype]["conv_ext"]],
#                                allow_missing=True)
#
#    return read_filelist

def get_main_read_filelist_for_correction(wildcards):
    read_filelist = []
    for main_datatype in assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options']["main_datatypes"]:
        if main_datatype not in config["data"]:
            continue
        read_filelist += expand(config["out_dir"]/ "data/{datatype}/final/{fileprefix}{extension}",
                                fileprefix=config["data"][main_datatype]["conv_file_prefix_list"],
                                extension=[config["data"][main_datatype]["conv_ext"]],
                                datatype=[main_datatype,],
                                allow_missing=True)
    return read_filelist

def get_read_filelist(wildcards, option_set_entry):
    read_filelist = []
    for datatype in stage_dict["contig"].parameters[wildcards.parameters]["option_set"][option_set_entry]:
        if datatype in config["data"]:
            read_filelist += expand(config["out_dir"]/ "data/{datatype}/final/{fileprefix}{extension}",
                                    fileprefix=config["data"][datatype]["conv_file_prefix_list"],
                                    datatype=[datatype,],
                                    extension=[config["data"][datatype]["conv_ext"]],
                                    allow_missing=True)
    return read_filelist

get_main_read_filelist = partial(get_read_filelist, option_set_entry="main_datatypes")
get_ultralong_read_filelist = partial(get_read_filelist, option_set_entry="ultra_long_reads")
get_hifi_read_filelist = partial(get_read_filelist, option_set_entry="hifi_datatypes")
get_nano_read_filelist = partial(get_read_filelist, option_set_entry="nano_datatypes")

#def get_ultralong_read_files(option_set):
#    read_filelist = []
#    for ultralong_read_type in option_set["ultra_long_reads"]:
#        if ultralong_read_type in config["data"]:
#            read_filelist += expand(config["out_dir"] / "data/{datatype}/final/{fileprefix}{extension}",
#                                    datatype=[ultralong_read_type,],
#                                    fileprefix=config["data"][ultralong_read_type]["conv_file_prefix_list"],
#                                    extension=[config["data"][ultralong_read_type]["conv_ext"]])
#    return read_filelist