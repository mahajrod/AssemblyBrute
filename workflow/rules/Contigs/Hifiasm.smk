#ruleorder: hifiasm_hic > hifiasm_hifi
ruleorder: hifiasm_hic > hifiasm_hic_tetra
localrules: get_lowcoverage_contig_ids, extract_lambda_value

def get_main_read_filelist_for_correction(wildcards):
    read_filelist = []
    for datatype in assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options']["main_datatypes"]:
        if datatype not in input_filedict:
            continue
        read_filelist += expand(output_dict["data"] / ("fastq/{datatype}/filtered/{fileprefix}%s" % config["fastq_extension"]),
                                fileprefix=input_file_prefix_dict[datatype],
                                datatype=[datatype,],
                                allow_missing=True)
    return read_filelist

rule hifiasm_correct:
    priority: 2000
    input:
        main_reads=get_main_read_filelist_for_correction,
        #hifi=expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #            fileprefix=input_file_prefix_dict["hifi"],
        #            allow_missing=True),
        #nanopore=expand(output_dict["data"] / ("fastq/nanopore/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                fileprefix=input_file_prefix_dict["nanopore"],
        #                allow_missing=True) if "nanopore" in input_filedict else [],
        #duplex=expand(output_dict["data"] / ("fastq/duplex/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                fileprefix=input_file_prefix_dict["duplex"],
        #                allow_missing=True) if "duplex" in input_filedict else [],
        #simplex=expand(output_dict["data"] / ("fastq/simplex/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                fileprefix=input_file_prefix_dict["simplex"],
        #                allow_missing=True) if "simplex" in input_filedict else [],
    output:
        ec_bin=output_dict["error_correction"] / "hifiasm_{correction_options, [^/]+}/{genome_prefix, [^/]+}.contig.ec.bin",
        ec_fasta=output_dict["error_correction"] / "hifiasm_{correction_options, [^/]+}/{genome_prefix, [^/]+}.contig.ec.fasta.gz",
        alias_ec_fasta=out_dir_path / "data/fastq/hifi/error_corrected_hifiasm_{correction_options, [^/]+}/{genome_prefix, [^/]+}.contig.ec.fasta.gz",
        ovlp_reverse_bin=output_dict["error_correction"] / "hifiasm_{correction_options, [^/]+}/{genome_prefix, [^/]+}.contig.ovlp.reverse.bin",
        ovlp_source_bin=output_dict["error_correction"] / "hifiasm_{correction_options, [^/]+}/{genome_prefix, [^/]+}.contig.ovlp.source.bin",
    params:
        window_size=lambda wildcards: parse_option("window_size", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -w "),
        bloom_filter_bits=lambda wildcards: parse_option("bloom_filter_bits", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -f "),
        rounds_of_error_correction=lambda wildcards: parse_option("rounds_of_error_correction", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -r "),
        length_of_adapters=lambda wildcards: parse_option("length_of_adapters", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -z "),
        max_kocc=lambda wildcards: parse_option("max-kocc", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " --max-kocc "),
        hg_size=lambda wildcards: parse_option("hg-size", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " --hg-size "),
        kmer_length=lambda wildcards: parse_option("kmer_len", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -k "),
        D=lambda wildcards: parse_option("D", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -D "), #" -D {0} ".format(parameters["tool_options"]["hifiasm"][wildcards.contig_options]["D"]) if "D" in parameters["tool_options"]["hifiasm"][wildcards.contig_options] else "",
        N=lambda wildcards: parse_option("N", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " -N "), #" -N {0} ".format(parameters["tool_options"]["hifiasm"][wildcards.contig_options]["N"]) if "N" in parameters["tool_options"]["hifiasm"][wildcards.contig_options] else "",
        ont_assembly=lambda wildcards: parse_option_flag("ont_mode", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " --ont "),
        #telomere_motif=lambda wildcards: parse_option("telomere_motif", config, " --telo-m ")
        #nanopore=(" --ul " + ",".join(map(str, expand(output_dict["data"] / ("fastq/nanopore/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                                              fileprefix=input_file_prefix_dict["nanopore"],
        #                                              allow_missing=True)))) if "nanopore" in input_filedict else "",
        #ultralong_reads=lambda wildcards: get_ultralong_read_files(input_file_prefix_dict,
        #                                                           stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set"]),
        #ul_cut=lambda wildcards: parse_option("ul-cut", assembler_option_set_group_dict["hifiasm"][wildcards.correction_options]['grouping_options'], " --ul-cut ")
    log:
        std=output_dict["log"] / "hifiasm_correct.{correction_options}.{genome_prefix}.log",
        pigz=output_dict["log"] / "hifiasm_correct.{correction_options}.{genome_prefix}.pigz.log",
        mv=output_dict["log"] / "hifiasm_correct.{correction_options}.{genome_prefix}.mv.log",
        ln=output_dict["log"] / "hifiasm_correct.{correction_options}.{genome_prefix}.ln.log",
        cluster_log=output_dict["cluster_log"] / "hifiasm_correct{correction_options}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "hifiasm_correct{correction_options}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "hifiasm_correct.{correction_options}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["hifiasm"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["hifiasm"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("hifiasm_correct"),
        cpus=parameters["threads"]["hifiasm"],
        time=parameters["time"]["hifiasm"],
        mem=parameters["memory_mb"]["hifiasm"],
    threads:
        parameters["threads"]["hifiasm"]
    shell:
         " OUTPUT_PREFIX={output.ec_bin}; "
         " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.ec.bin}}; "
         " UNCOMPRESSED_FASTA={output.ec_fasta}; "
         " UNCOMPRESSED_FASTA=${{UNCOMPRESSED_FASTA%sta.gz}}; "
         " hifiasm -t {threads} -e --write-ec {params.window_size} {params.bloom_filter_bits} {params.ont_assembly} "
         " {params.rounds_of_error_correction} {params.length_of_adapters} {params.max_kocc} {params.hg_size}"
         " {params.kmer_length} {params.D} {params.N} "
         " -o ${{OUTPUT_PREFIX}} {input.main_reads}  1>{log.std} 2>&1;"
         " pigz -p {threads} ${{UNCOMPRESSED_FASTA}} > {log.pigz} 2>&1 ; "
         " mv ${{UNCOMPRESSED_FASTA}}.gz {output.ec_fasta} > {log.mv} 2>&1; "
         " ln -sf ../../../../../{output.ec_fasta} {output.alias_ec_fasta} > {log.ln} 2>&1; "
         " sleep 60; "

def get_ultralong_read_files(input_file_prefix_dict, option_set):
    #option_set = stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set"]
    read_filelist = []
    for ultralong_read_type in option_set["ultra_long_reads"]:
        if ultralong_read_type in input_file_prefix_dict:
            read_filelist += expand(output_dict["data"] / ("fastq/%s/filtered/{fileprefix}%s" % (ultralong_read_type,
                                                                                                 config["fastq_extension"])),
                                    fileprefix=input_file_prefix_dict[ultralong_read_type])
    return read_filelist

def get_coverage_estimator(wildcards):
    coverage_estimator = parameters["tool_options"]["hifiasm"][wildcards.contig_options]["coverage_estimator"]
    if "hifiasm" in config["tool_manually_adjusted_features"]:
        if "coverage_estimator" in config["tool_manually_adjusted_features"]["hifiasm"]:
            if config["tool_manually_adjusted_features"]["hifiasm"]["coverage_estimator"] is not None:
                coverage_estimator = config["tool_manually_adjusted_features"]["hifiasm"]["coverage_estimator"]

    if coverage_estimator not in config["allowed_coverage_estimators"]:
        raise ValueError("ERROR!!! Unrecognized coverage estimator ({0}) for contig assembly hifiasm_{1} !".format(coverage_estimator,
                                                                                                                   wildcards.contig_options))
    return coverage_estimator

def get_coverage_estimator_report_filename(wildcards):
    coverage_estimator = get_coverage_estimator(wildcards)
    report_filename = output_dict["kmer"] / ("{0}/filtered/{1}/{2}.{3}.filtered.{4}.{5}.{1}.parameters".format("_".join(parameters["tool_options"]["hifiasm"][wildcards.contig_options]["main_datatypes"]), #config["final_kmer_datatype"]
                                                                                                               coverage_estimator,
                                                                                                               wildcards.genome_prefix,
                                                                                                               "_".join(parameters["tool_options"]["hifiasm"][wildcards.contig_options]["main_datatypes"]), #config["final_kmer_datatype"],
                                                                                                               config["final_kmer_length"],
                                                                                                               config["final_kmer_counter"]))
    #print(report_filename)
    #print(coverage_estimator)
    return report_filename

"""
def get_lambda_value(wildcards):
    coverage_estimator = get_coverage_estimator(wildcards)
    if "hifiasm" in contig["tool_manually_adjusted_features"]:
        if "lambda" in contig["tool_manually_adjusted_features"]["hifiasm"]:
            if contig["tool_manually_adjusted_features"]["hifiasm"]["lambda"] is not None:
                if isinstance(contig["tool_manually_adjusted_features"]["hifiasm"]["lambda"], Number):
                    lambda_value = contig["tool_manually_adjusted_features"]["hifiasm"]["lambda"]
                    print("Using a preset lambda value ({0}) for contig assembly hifiasm_{1} ...".format(contig["tool_manually_adjusted_features"]["hifiasm"]["lambda"],
                                                                                                         wildcards.contig_options))
                    return lambda_value
                else:
                    raise ValueError("ERROR!!! Preset lambda value is not a number! Check value in contig['tool_manually_adjusted_features']['hifiasm']['lambda'] ...")

    print("Using the {0} as a coverage estimator for contig assembly hifiasm_{1} ...".format(coverage_estimator,
                                                                                           wildcards.contig_options))
    report_filename = output_dict["kmer"] / ("{0}/filtered/{1}/{2}.{3}.filtered.{4}.{5}.{1}.parameters".format(config["final_kmer_datatype"],
                                                                                                               coverage_estimator,
                                                                                                               wildcards.genome_prefix,
                                                                                                               config["final_kmer_datatype"],
                                                                                                               config["final_kmer_length"],
                                                                                                               config["final_kmer_counter"]))
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
        raise ValueError("ERROR!!! Unimplemented (but planned) coverage estimator ({0}) for contig assembly hifiasm_{1} !".format(coverage_estimator,
                                                                                                                                  wildcards.contig_options))
    return lambda_value
"""
rule extract_lambda_value:
    input:
        coverage_estimator_report_filename=get_coverage_estimator_report_filename
    output:
        lambda_file=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.lambda",
    log:
        std=output_dict["log"] / "extract_lambda_value.{contig_options}.{genome_prefix}.log",
        cluster_log=output_dict["cluster_log"] / "extract_lambda_value.{contig_options}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "extract_lambda_value.{contig_options}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "extract_lambda_value.{contig_options}.{genome_prefix}.benchmark.txt"
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("extract_lambda_value"),
        cpus=parameters["threads"]["extract_lambda_value"],
        time=parameters["time"]["extract_lambda_value"],
        mem=parameters["memory_mb"]["extract_lambda_value"],
    threads:
        parameters["threads"]["extract_lambda_value"]
    run:
        with open(log.std, "w") as log_fd, open(output.lambda_file, "w") as out_fd:
            coverage_estimator = get_coverage_estimator(wildcards)
            if "hifiasm" in config["tool_manually_adjusted_features"]:
                if "lambda" in config["tool_manually_adjusted_features"]["hifiasm"]:
                    if config["tool_manually_adjusted_features"]["hifiasm"]["lambda"] is not None:
                        if isinstance(config["tool_manually_adjusted_features"]["hifiasm"]["lambda"], Number):
                            lambda_value = config["tool_manually_adjusted_features"]["hifiasm"]["lambda"]
                            print("Using a preset lambda value ({0}) for contig assembly hifiasm_{1} ...".format(config["tool_manually_adjusted_features"]["hifiasm"]["lambda"],
                                                                                                                 wildcards.contig_options))
                            return lambda_value
                        else:
                            message = "ERROR!!! Preset lambda value is not a number! Check value in contig['tool_manually_adjusted_features']['hifiasm']['lambda'] ..."
                            log_fd.write(message + "\n")
                            raise ValueError(message)

            log_fd.write("Using the {0} as a coverage estimator for contig assembly hifiasm_{1} ...\n".format(coverage_estimator,
                                                                                                            wildcards.contig_options))
            report_filename = output_dict["kmer"] / ("{0}/filtered/{1}/{2}.{3}.filtered.{4}.{5}.{1}.parameters".format("_".join(parameters["tool_options"]["hifiasm"][wildcards.contig_options]["main_datatypes"]), #config["final_kmer_datatype"],
                                                                                                                       coverage_estimator,
                                                                                                                       wildcards.genome_prefix,
                                                                                                                       "_".join(parameters["tool_options"]["hifiasm"][wildcards.contig_options]["main_datatypes"]), #config["final_kmer_datatype"],
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
                                                                                                                                   wildcards.contig_options)
                log_fd.write(message + "\n")
                raise ValueError(message)
            log_fd.write("Lambda:\t%.2f\n" % lambda_value)
            out_fd.write("%.2f\n" % lambda_value)

def get_main_read_filelist(wildcards):
    read_filelist = []
    for datatype in parameters["tool_options"]["hifiasm"][wildcards.contig_options]["main_datatypes"]:
        if datatype not in input_filedict:
            continue
        read_filelist += expand(output_dict["data"] / ("fastq/{datatype}/filtered/{fileprefix}%s" % config["fastq_extension"]),
                                fileprefix=input_file_prefix_dict[datatype],
                                datatype=[datatype,],
                                allow_missing=True)
    #print("AAAA")
    #print(wildcards)
    #print(read_filelist)
    return read_filelist

rule hifiasm_hic: # TODO: add support for polyploid assemblies
    priority: 1000
    input:
        main_reads=get_main_read_filelist,
        #hifi=expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #            fileprefix=input_file_prefix_dict["hifi"],
        #            allow_missing=True),
        ultralong_reads=lambda wildcards: get_ultralong_read_files(input_file_prefix_dict,
                                                                   stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set"]),
        #nanopore=expand(output_dict["data"] / ("fastq/nanopore/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                fileprefix=input_file_prefix_dict["nanopore"],
        #                allow_missing=True) if "nanopore" in input_filedict else [],
        #lqccs=expand(output_dict["data"] / ("fastq/lqccs/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                fileprefix=input_file_prefix_dict["lqccs"],
        #                allow_missing=True) if "lqccs" in input_filedict else [],
        #hic_forward=list(map(lambda s: output_dict["data"] / "fastq/hic/raw/" / s.name, input_filedict["hic"][::2])) if "hic" in input_filedict else [],
        #hic_reverse=list(map(lambda s: output_dict["data"] / "fastq/hic/raw/" / s.name, input_filedict["hic"][1::2])) if "hic" in input_filedict else [],
        hic_forward=expand(output_dict["data"] / ("fastq/hic/filtered/{pairprefix}_1%s" % config["fastq_extension"]), pairprefix=input_pairprefix_dict["hic"]) if "hic" in input_filedict else [],
        hic_reverse=expand(output_dict["data"] / ("fastq/hic/filtered/{pairprefix}_2%s" % config["fastq_extension"]), pairprefix=input_pairprefix_dict["hic"]) if "hic" in input_filedict else [],
        ec_bin=lambda wildcards: output_dict["error_correction"] / "hifiasm_{0}/{1}.contig.ec.bin".format(stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set_group"],
                                                                                                          wildcards.genome_prefix),
        ovlp_reverse_bin=lambda wildcards: output_dict["error_correction"] / "hifiasm_{0}/{1}.contig.ovlp.reverse.bin".format(stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set_group"],
                                                                                                                              wildcards.genome_prefix),
        ovlp_source_bin=lambda wildcards: output_dict["error_correction"] / "hifiasm_{0}/{1}.contig.ovlp.source.bin".format(stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set_group"],
                                                                                                                            wildcards.genome_prefix),
        #coverage_estimator_report_filename=get_coverage_estimator_report_filename
        lambda_file=rules.extract_lambda_value.output.lambda_file
    output:
        primary_contig_graph=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hic.hap1.p_ctg.gfa",
        alternative_contig_graph=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hic.hap2.p_ctg.gfa",
        alt_contig_graph=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hic.a_ctg.gfa",
        primary_alias=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hap1.unfiltered.gfa",
        alternative_alias=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hap2.unfiltered.gfa",
        alt_alias=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.alt.unfiltered.gfa",
    params:
        purge_level=lambda wildcards: parameters["tool_options"]["hifiasm"][wildcards.contig_options]["purge level"],
        ploidy=lambda wildcards: stage_dict["contig"]["parameters"][f"hifiasm_{wildcards.contig_options}"]["option_set"]["assembly_ploidy"], #config["ploidy"],
        cov_multiplicator=lambda wildcards: parameters["tool_options"]["hifiasm"][wildcards.contig_options]["cov_multiplicator"],
        window_size=lambda wildcards: parse_option("window_size", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -w "),
        bloom_filter_bits=lambda wildcards: parse_option("bloom_filter_bits", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -f "),
        rounds_of_error_correction=lambda wildcards: parse_option("rounds_of_error_correction", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -r "),
        length_of_adapters=lambda wildcards: parse_option("length_of_adapters", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -z "),
        max_kocc=lambda wildcards: parse_option("max-kocc", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --max-kocc "),
        hg_size=lambda wildcards: parse_option("hg-size", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --hg-size "),
        kmer_length=lambda wildcards: parse_option("kmer_len", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -k "),
        dual_scaf=lambda wildcards: parse_option_flag("dual_scaf", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --dual-scaf "),
        D=lambda wildcards: parse_option("D", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -D "), #" -D {0} ".format(parameters["tool_options"]["hifiasm"][wildcards.contig_options]["D"]) if "D" in parameters["tool_options"]["hifiasm"][wildcards.contig_options] else "",
        N=lambda wildcards: parse_option("N", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -N "),
        ignore_bin=lambda wildcards: " -i " if ("ignore_bin" in parameters["tool_options"]["hifiasm"][wildcards.contig_options]) and parameters["tool_options"]["hifiasm"][wildcards.contig_options]["ignore_bin"] else "",
        hic_forward=(" --h1 " + ",".join(map(str, expand(output_dict["data"] / ("fastq/hic/filtered/{pairprefix}_1%s" % config["fastq_extension"]), pairprefix=input_pairprefix_dict["hic"]) ))) if "hic" in input_filedict else "", #in case of multiple hic libraries files in the list MUST be COMMA-separated
        hic_reverse=(" --h2 " + ",".join(map(str, expand(output_dict["data"] / ("fastq/hic/filtered/{pairprefix}_2%s" % config["fastq_extension"]), pairprefix=input_pairprefix_dict["hic"]) ))) if "hic" in input_filedict else "",
        ultralong_reads=lambda wildcards: (" --ul " + ",".join(map(str,
                                                                   get_ultralong_read_files(input_file_prefix_dict,
                                                                                            stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set"])
                                                                   )
                                                              )
                                           ) if get_ultralong_read_files(input_file_prefix_dict,
                                                                         stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set"]) else "",
        telomere_motif= lambda wildcards: parse_option("telomere_motif", config, " --telo-m ") if parameters["tool_options"]["hifiasm"][wildcards.contig_options]["use_telomere"] else "",
        #nanopore=(" --ul " + ",".join(map(str, expand(output_dict["data"] / ("fastq/nanopore/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                                              fileprefix=input_file_prefix_dict["nanopore"],
        #                                              allow_missing=True)))) if "nanopore" in input_filedict else "",
        #lqccs=(" --ul " + ",".join(map(str, expand(output_dict["data"] / ("fastq/lqccs/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                                              fileprefix=input_file_prefix_dict["nanopore"],
        #                                              allow_missing=True)))) if "nanopore" in input_filedict else "",
        ul_cut=lambda wildcards: parse_option("ul-cut", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --ul-cut "),
        ont_assembly= lambda wildcards: parse_option_flag("ont_mode", parameters["tool_options"]["hifiasm"][wildcards.contig_options]," --ont "),
    log:
        std=output_dict["log"] / "hifiasm.{contig_options}.{genome_prefix}.log",
        cluster_log=output_dict["cluster_log"] / "hifiasm.{contig_options}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "hifiasm.{contig_options}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "hifiasm.{contig_options}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["hifiasm"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["hifiasm"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("hifiasm_hic"),
        cpus=parameters["threads"]["hifiasm"],
        time=parameters["time"]["hifiasm"],
        mem=parameters["memory_mb"]["hifiasm"],
    threads:
        parameters["threads"]["hifiasm"]
    shell:
         " OUTPUT_PREFIX={output.primary_alias}; "
         " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.hap1.unfiltered.gfa}}; "
         " OUT_DIR=`dirname ${{OUTPUT_PREFIX}}`; "
         " ln -sf ../../../{input.ec_bin} ${{OUT_DIR}} 1>{log.std} 2>&1; "
         " ln -sf ../../../{input.ovlp_reverse_bin} ${{OUT_DIR}} 1>>{log.std} 2>&1; "
         " ln -sf ../../../{input.ovlp_source_bin} ${{OUT_DIR}} 1>>{log.std} 2>&1; "
         " LAMBDA=`head -n 1 {input.lambda_file}` 1>>{log.std} 2>&1;  "
         " COV_UPPER_BOUNDARY=`echo \"{params.cov_multiplicator}*${{LAMBDA}}\" | bc` 1>>{log.std} 2>&1;  "
         " COV_UPPER_BOUNDARY=${{COV_UPPER_BOUNDARY%.*}}; "
         " hifiasm {params.window_size} {params.bloom_filter_bits} {params.ont_assembly} "
         " {params.rounds_of_error_correction} {params.length_of_adapters} {params.max_kocc} {params.hg_size}"
         " {params.kmer_length} {params.D} {params.N} {params.ignore_bin} --primary -t {threads} -l {params.purge_level}  -o ${{OUTPUT_PREFIX}} "
         " --n-hap {params.ploidy} --purge-max ${{COV_UPPER_BOUNDARY}} "
         " {params.hic_forward} {params.hic_reverse} {params.ultralong_reads} {params.ul_cut} {params.dual_scaf} "
         " {params.telomere_motif} "
         " {input.main_reads}  1>{log.std} 2>&1; "         
         " ln -sf `basename {output.primary_contig_graph}` {output.primary_alias} 1>>{log.std} 2>&1; "
         " ln -sf `basename {output.alternative_contig_graph}` {output.alternative_alias} 1>>{log.std} 2>&1; "
         " ln -sf `basename {output.alt_contig_graph}` {output.alt_alias} 1>>{log.std} 2>&1; "
         " sleep 60;"
         #" COV_UPPER_BOUNDARY=`awk 'NR==2 {{printf \"%.0f\", {params.cov_multiplicator} * $2}}' {input.genomescope_report}`; "

rule hifiasm_hic_tetra: # TODO: add support for polyploid assemblies
    priority: 1000
    input:
        main_reads=get_main_read_filelist,
        #hifi=expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #            fileprefix=input_file_prefix_dict["hifi"],
        #            allow_missing=True),
        ultralong_reads=lambda wildcards: get_ultralong_read_files(input_file_prefix_dict,
                                                                   stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set"]),
        #nanopore=expand(output_dict["data"] / ("fastq/nanopore/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                fileprefix=input_file_prefix_dict["nanopore"],
        #                allow_missing=True) if "nanopore" in input_filedict else [],
        #lqccs=expand(output_dict["data"] / ("fastq/lqccs/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                fileprefix=input_file_prefix_dict["lqccs"],
        #                allow_missing=True) if "lqccs" in input_filedict else [],
        #hic_forward=list(map(lambda s: output_dict["data"] / "fastq/hic/raw/" / s.name, input_filedict["hic"][::2])) if "hic" in input_filedict else [],
        #hic_reverse=list(map(lambda s: output_dict["data"] / "fastq/hic/raw/" / s.name, input_filedict["hic"][1::2])) if "hic" in input_filedict else [],
        hic_forward=expand(output_dict["data"] / ("fastq/hic/filtered/{pairprefix}_1%s" % config["fastq_extension"]), pairprefix=input_pairprefix_dict["hic"]) if "hic" in input_filedict else [],
        hic_reverse=expand(output_dict["data"] / ("fastq/hic/filtered/{pairprefix}_2%s" % config["fastq_extension"]), pairprefix=input_pairprefix_dict["hic"]) if "hic" in input_filedict else [],
        ec_bin=lambda wildcards: output_dict["error_correction"] / "hifiasm_{0}/{1}.contig.ec.bin".format(stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set_group"],
                                                                                                          wildcards.genome_prefix),
        ovlp_reverse_bin=lambda wildcards: output_dict["error_correction"] / "hifiasm_{0}/{1}.contig.ovlp.reverse.bin".format(stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set_group"],
                                                                                                                              wildcards.genome_prefix),
        ovlp_source_bin=lambda wildcards: output_dict["error_correction"] / "hifiasm_{0}/{1}.contig.ovlp.source.bin".format(stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set_group"],
                                                                                                                            wildcards.genome_prefix),
        #coverage_estimator_report_filename=get_coverage_estimator_report_filename
        lambda_file=rules.extract_lambda_value.output.lambda_file
    output:
        primary_contig_graph=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hic.hap1.p_ctg.gfa",
        alternative_contig_graph=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hic.hap2.p_ctg.gfa",
        third_contig_graph=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hic.hap3.p_ctg.gfa",
        forth_contig_graph=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hic.hap4.p_ctg.gfa",
        #alt_contig_graph=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hic.a_ctg.gfa",
        primary_alias=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hap1.unfiltered.gfa",
        alternative_alias=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hap2.unfiltered.gfa",
        third_alias=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hap3.unfiltered.gfa",
        forth_alias=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hap4.unfiltered.gfa",
        #alt_alias=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.alt.unfiltered.gfa",

        #contig_graph_list=expand(output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.hic.{haplotype}.p_ctg.gfa",
        #                         haplotype=["hap1", "hap2", "hap3", "hap4",], allow_missing=True),
        alt_graph=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hic.a_ctg.gfa",
        #contig_graph_alias_list=expand(output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.{haplotype}.unfiltered.gfa",
        #                         haplotype=["hap1", "hap2", "hap3", "hap4"], allow_missing=True),
        alt_alias=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.alt.unfiltered.gfa",
    params:
        purge_level=lambda wildcards: parameters["tool_options"]["hifiasm"][wildcards.contig_options]["purge level"],
        ploidy=lambda wildcards: stage_dict["contig"]["parameters"][f"hifiasm_{wildcards.contig_options}"]["option_set"]["assembly_ploidy"], #config["ploidy"],
        cov_multiplicator=lambda wildcards: parameters["tool_options"]["hifiasm"][wildcards.contig_options]["cov_multiplicator"],
        window_size=lambda wildcards: parse_option("window_size", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -w "),
        bloom_filter_bits=lambda wildcards: parse_option("bloom_filter_bits", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -f "),
        rounds_of_error_correction=lambda wildcards: parse_option("rounds_of_error_correction", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -r "),
        length_of_adapters=lambda wildcards: parse_option("length_of_adapters", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -z "),
        max_kocc=lambda wildcards: parse_option("max-kocc", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --max-kocc "),
        hg_size=lambda wildcards: parse_option("hg-size", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --hg-size "),
        kmer_length=lambda wildcards: parse_option("kmer_len", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -k "),
        dual_scaf=lambda wildcards: parse_option_flag("dual_scaf", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --dual-scaf "),
        D=lambda wildcards: parse_option("D", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -D "), #" -D {0} ".format(parameters["tool_options"]["hifiasm"][wildcards.contig_options]["D"]) if "D" in parameters["tool_options"]["hifiasm"][wildcards.contig_options] else "",
        N=lambda wildcards: parse_option("N", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -N "),
        ignore_bin=lambda wildcards: " -i " if ("ignore_bin" in parameters["tool_options"]["hifiasm"][wildcards.contig_options]) and parameters["tool_options"]["hifiasm"][wildcards.contig_options]["ignore_bin"] else "",
        hic_forward=(" --h1 " + ",".join(map(str, expand(output_dict["data"] / ("fastq/hic/filtered/{pairprefix}_1%s" % config["fastq_extension"]), pairprefix=input_pairprefix_dict["hic"]) ))) if "hic" in input_filedict else "", #in case of multiple hic libraries files in the list MUST be COMMA-separated
        hic_reverse=(" --h2 " + ",".join(map(str, expand(output_dict["data"] / ("fastq/hic/filtered/{pairprefix}_2%s" % config["fastq_extension"]), pairprefix=input_pairprefix_dict["hic"]) ))) if "hic" in input_filedict else "",
        ultralong_reads=lambda wildcards: (" --ul " + ",".join(map(str,
                                                                   get_ultralong_read_files(input_file_prefix_dict,
                                                                                            stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set"])
                                                                   )
                                                              )
                                           ) if get_ultralong_read_files(input_file_prefix_dict,
                                                                         stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set"]) else "",
        telomere_motif= lambda wildcards: parse_option("telomere_motif", config, " --telo-m ") if parameters["tool_options"]["hifiasm"][wildcards.contig_options]["use_telomere"] else "",
        #nanopore=(" --ul " + ",".join(map(str, expand(output_dict["data"] / ("fastq/nanopore/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                                              fileprefix=input_file_prefix_dict["nanopore"],
        #                                              allow_missing=True)))) if "nanopore" in input_filedict else "",
        #lqccs=(" --ul " + ",".join(map(str, expand(output_dict["data"] / ("fastq/lqccs/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                                              fileprefix=input_file_prefix_dict["nanopore"],
        #                                              allow_missing=True)))) if "nanopore" in input_filedict else "",
        ul_cut=lambda wildcards: parse_option("ul-cut", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --ul-cut "),
        ont_assembly= lambda wildcards: parse_option_flag("ont_mode", parameters["tool_options"]["hifiasm"][wildcards.contig_options]," --ont "),
    log:
        std=output_dict["log"] / "hifiasm.{contig_options}.{genome_prefix}.log",
        cluster_log=output_dict["cluster_log"] / "hifiasm.{contig_options}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "hifiasm.{contig_options}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "hifiasm.{contig_options}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["hifiasm"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["hifiasm"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("hifiasm_hic"),
        cpus=parameters["threads"]["hifiasm"],
        time=parameters["time"]["hifiasm"],
        mem=parameters["memory_mb"]["hifiasm"],
    threads:
        parameters["threads"]["hifiasm"]
    shell:
         " OUTPUT_PREFIX={output.alt_alias}; "
         " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.alt.unfiltered.gfa}}; "
         " OUT_DIR=`dirname ${{OUTPUT_PREFIX}}`; "
         " ln -sf ../../../{input.ec_bin} ${{OUT_DIR}} 1>{log.std} 2>&1; "
         " ln -sf ../../../{input.ovlp_reverse_bin} ${{OUT_DIR}} 1>>{log.std} 2>&1; "
         " ln -sf ../../../{input.ovlp_source_bin} ${{OUT_DIR}} 1>>{log.std} 2>&1; "
         " LAMBDA=`head -n 1 {input.lambda_file}` 1>>{log.std} 2>&1;  "
         " COV_UPPER_BOUNDARY=`echo \"{params.cov_multiplicator}*${{LAMBDA}}\" | bc` 1>>{log.std} 2>&1;  "
         " COV_UPPER_BOUNDARY=${{COV_UPPER_BOUNDARY%.*}}; "
         " hifiasm {params.window_size} {params.bloom_filter_bits} {params.ont_assembly} "
         " {params.rounds_of_error_correction} {params.length_of_adapters} {params.max_kocc} {params.hg_size}"
         " {params.kmer_length} {params.D} {params.N} {params.ignore_bin} --primary -t {threads} -l {params.purge_level}  -o ${{OUTPUT_PREFIX}} "
         " --n-hap {params.ploidy} --purge-max ${{COV_UPPER_BOUNDARY}} "
         " {params.hic_forward} {params.hic_reverse} {params.ultralong_reads} {params.ul_cut} {params.dual_scaf} "
         " {params.telomere_motif} "
         " {input.main_reads}  1>{log.std} 2>&1; "         
         " ln -sf `basename {output.alt_graph}` {output.alt_alias} 1>>{log.std} 2>&1; "
         " OUTPUT_PREFIX_BASENAME=`basename ${{OUTPUT_PREFIX}}`; "
         " for HAP in hap1 hap2 hap3 hap4; "
         "    do "
         "    ln -sf ${{OUTPUT_PREFIX_BASENAME}}.hic.${{HAP}}.p_ctg.gfa ${{OUTPUT_PREFIX_BASENAME}}.${{HAP}}.unfiltered.gfa >>{log.std}  2>&1;  "
         "    done; "
         " sleep 60; "


rule hifiasm_long_reads_only:
    priority: 1000
    input:
        main_reads=get_main_read_filelist,
        #hifi=expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #            fileprefix=input_file_prefix_dict["hifi"],
        #            allow_missing=True),
        #nanopore=expand(output_dict["data"] / ("fastq/nanopore/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                fileprefix=input_file_prefix_dict["nanopore"],
        #                allow_missing=True) if "nanopore" in input_filedict else [],
        #lqccs=expand(output_dict["data"] / ("fastq/lqccs/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                fileprefix=input_file_prefix_dict["lqccs"],
        #                allow_missing=True) if "lqccs" in input_filedict else [],
        ultralong_reads=lambda wildcards: get_ultralong_read_files(input_file_prefix_dict,
                                                                   stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set"]),
        ec_bin=lambda wildcards: output_dict["error_correction"] / "hifiasm_{0}/{1}.contig.ec.bin".format(stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set_group"],
                                                                                                          wildcards.genome_prefix),
        ovlp_reverse_bin=lambda wildcards: output_dict["error_correction"] / "hifiasm_{0}/{1}.contig.ovlp.reverse.bin".format(stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set_group"],
                                                                                                                              wildcards.genome_prefix),
        ovlp_source_bin=lambda wildcards: output_dict["error_correction"] / "hifiasm_{0}/{1}.contig.ovlp.source.bin".format(stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set_group"],
                                                                                                                            wildcards.genome_prefix),
        lambda_file=rules.extract_lambda_value.output.lambda_file
    output:
        primary_contig_graph=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.p_ctg.gfa",
        alt_contig_graph=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.a_ctg.gfa",
        primary_alias=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.hap0.unfiltered.gfa",
        alt_alias=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.alt0.unfiltered.gfa",
    params:
        purge_level=lambda wildcards: parameters["tool_options"]["hifiasm"][wildcards.contig_options]["purge level"],
        ploidy=lambda wildcards: stage_dict["contig"]["parameters"][f"hifiasm_{wildcards.contig_options}"]["option_set"]["assembly_ploidy"], #config["ploidy"],
        cov_multiplicator=lambda wildcards: parameters["tool_options"]["hifiasm"][wildcards.contig_options]["cov_multiplicator"],
        window_size=lambda wildcards: parse_option("window_size", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -w "),
        bloom_filter_bits=lambda wildcards: parse_option("bloom_filter_bits", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -f "),
        rounds_of_error_correction=lambda wildcards: parse_option("rounds_of_error_correction", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -r "),
        length_of_adapters=lambda wildcards: parse_option("length_of_adapters", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -z "),
        max_kocc=lambda wildcards: parse_option("max-kocc", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --max-kocc "),
        hg_size=lambda wildcards: parse_option("hg-size", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --hg-size "),
        kmer_length=lambda wildcards: parse_option("kmer_len", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -k "),
        D=lambda wildcards: parse_option("D", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -D "), #" -D {0} ".format(parameters["tool_options"]["hifiasm"][wildcards.contig_options]["D"]) if "D" in parameters["tool_options"]["hifiasm"][wildcards.contig_options] else "",
        N=lambda wildcards: parse_option("N", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " -N "),
        dual_scaf=lambda wildcards: parse_option_flag("dual_scaf", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --dual-scaf "),
        ignore_bin=lambda wildcards: " -i " if ("ignore_bin" in parameters["tool_options"]["hifiasm"][wildcards.contig_options]) and parameters["tool_options"]["hifiasm"][wildcards.contig_options]["ignore_bin"] else "",
        telomere_motif= lambda wildcards: parse_option("telomere_motif",config," --telo-m ")  if parameters["tool_options"]["hifiasm"][wildcards.contig_options]["use_telomere"] else "",
        #nanopore=(" --ul " + ",".join(map(str, expand(output_dict["data"] / ("fastq/nanopore/filtered/{fileprefix}%s" % config["fastq_extension"]),
        #                                              fileprefix=input_file_prefix_dict["nanopore"],
        #                                              allow_missing=True)))) if "nanopore" in input_filedict else "",
        ultralong_reads=lambda wildcards: (" --ul " + ",".join(map(str,
                                                                   get_ultralong_read_files(input_file_prefix_dict,
                                                                                            stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set"])
                                                                   )
                                                              )
                                           ) if get_ultralong_read_files(input_file_prefix_dict,
                                                                         stage_dict["contig"]["parameters"]["hifiasm_" + wildcards.contig_options]["option_set"]) else "",
        ul_cut=lambda wildcards: parse_option("ul-cut", parameters["tool_options"]["hifiasm"][wildcards.contig_options], " --ul-cut "),
        ont_assembly= lambda wildcards: parse_option_flag("ont_mode", parameters["tool_options"]["hifiasm"][wildcards.contig_options]," --ont "),
    log:
        std=output_dict["log"] / "hifiasm.{contig_options}.{genome_prefix}.log",
        cluster_log=output_dict["cluster_log"] / "hifiasm.{contig_options}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "hifiasm.{contig_options}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "hifiasm.{contig_options}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["hifiasm"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["hifiasm"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("hifiasm_hifi"),
        cpus=parameters["threads"]["hifiasm"],
        time=parameters["time"]["hifiasm"],
        mem=parameters["memory_mb"]["hifiasm"],
    threads:
        parameters["threads"]["hifiasm"]
    shell:
         " OUTPUT_PREFIX={output.primary_alias}; "
         " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.hap0.unfiltered.gfa}}; "
         " OUT_DIR=`dirname ${{OUTPUT_PREFIX}}`; "
         " ln -sf ../../../{input.ec_bin} ${{OUT_DIR}}; "
         " ln -sf ../../../{input.ovlp_reverse_bin} ${{OUT_DIR}}; "
         " ln -sf ../../../{input.ovlp_source_bin} ${{OUT_DIR}}; "
         " LAMBDA=`head -n 1 {input.lambda_file}`; "
         " COV_UPPER_BOUNDARY=`echo \"{params.cov_multiplicator}*${{LAMBDA}}\" | bc`; "
         " COV_UPPER_BOUNDARY=${{COV_UPPER_BOUNDARY%.*}}; "
         " hifiasm {params.window_size} {params.bloom_filter_bits} {params.ont_assembly} "
         " {params.rounds_of_error_correction} {params.length_of_adapters} {params.max_kocc} {params.hg_size} "
         " {params.kmer_length} {params.D} {params.N} {params.ignore_bin} {params.ul_cut}"
         " --primary -t {threads} -l {params.purge_level}  -o ${{OUTPUT_PREFIX}} "
         " --n-hap {params.ploidy} --purge-max ${{COV_UPPER_BOUNDARY}} {params.ultralong_reads} {params.dual_scaf} "
         " {params.telomere_motif} "
         " {input.main_reads}  1>{log.std} 2>&1;"
         " ln -sf `basename {output.primary_contig_graph}` {output.primary_alias};"
         " ln -sf `basename {output.alt_contig_graph}` {output.alt_alias};"

         #" COV_UPPER_BOUNDARY=`awk 'NR==2 {{printf \"%.0f\", {params.cov_multiplicator} * $2}}' {input.genomescope_report}`; "

rule get_lowcoverage_contig_ids:
    input:
        cov=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.{haplotype}.unfiltered.gfa.cov"
    output:
        low_cov_ids=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.{haplotype, [^/]+}.unfiltered.gfa.lowcov.ids",
    log:
        std=output_dict["log"] / "get_lowcoverage_contig_ids.{contig_options}.{genome_prefix}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "get_lowcoverage_contig_ids.{contig_options}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "get_lowcoverage_contig_ids.{contig_options}.{genome_prefix}.{haplotype}.cluster.err",
    params:
        min_coverage=lambda wildcards: parameters["tool_options"]["hifiasm"][wildcards.contig_options]["min_contig_coverage"]
    benchmark:
        output_dict["benchmark"] / "get_lowcoverage_contig_ids.{contig_options}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("get_lowcoverage_contig_ids"),
        cpus=parameters["threads"]["get_lowcoverage_contig_ids"],
        time=parameters["time"]["get_lowcoverage_contig_ids"],
        mem=parameters["memory_mb"]["get_lowcoverage_contig_ids"],
    threads:
        parameters["threads"]["get_lowcoverage_contig_ids"]
    shell:
         " if [ '{wildcards.haplotype}' == 'alt' ];"
         " then "
         "   > {output.low_cov_ids}; "
         " else "
         "   awk '{{if ($2 < {params.min_coverage}) print $1}}' {input.cov} > {output.low_cov_ids} 2>{log.std}; "
         " fi; "

rule filter_contigs_by_coverage:
    input:
        low_cov_ids=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.{haplotype}.unfiltered.gfa.lowcov.ids",
        unfiltered_fasta=output_dict["contig"] / "hifiasm_{contig_options}/{genome_prefix}.contig.{haplotype}.unfiltered.fasta"
    output:
        filtered_fasta=output_dict["contig"] / "hifiasm_{contig_options, [^/]+}/{genome_prefix, [^/]+}.contig.{haplotype, [^/]+}.lenfiltered.fasta",
    log:
        std=output_dict["log"] / "filter_contigs_by_coverage.{contig_options}.{genome_prefix}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "filter_contigs_by_coverage.{contig_options}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "filter_contigs_by_coverage.{contig_options}.{genome_prefix}.{haplotype}.cluster.err",
    benchmark:
        output_dict["benchmark"] / "filter_contigs_by_coverage.{contig_options}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("filter_contigs_by_coverage"),
        cpus=parameters["threads"]["filter_contigs_by_coverage"],
        time=parameters["time"]["filter_contigs_by_coverage"],
        mem=parameters["memory_mb"]["filter_contigs_by_coverage"],
    threads:
        parameters["threads"]["filter_contigs_by_coverage"]
    shell:
         " extract_sequences_by_ids.py -i {input.unfiltered_fasta} -d {input.low_cov_ids} -r "
         " -o {output.filtered_fasta} > {log.std} 2>&1; "

