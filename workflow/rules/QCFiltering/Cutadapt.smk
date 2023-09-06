ruleorder: cutadapt_illumina > cutadapt > porechop_abi
rule cutadapt:
    input:
        fastq=output_dict["data"] / ("fastq/hifi/raw/{fileprefix}%s" % config["fastq_extension"])
    output:
        fastq=output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
        stats=output_dict["data"] / "fastq/hifi/filtered/{fileprefix}.cutadapt.stats"
    params:
        error_rate=lambda wildcards: "-e {0} ".format(parameters["tool_options"]["cutadapt"]["hifi"]["error_rate"]) if "error_rate" in parameters["tool_options"]["cutadapt"]["hifi"] else "",
        min_read_length=lambda wildcards: " -m {0} ".format(parameters["tool_options"]["cutadapt"]["hifi"]["min_read_length"]) if "min_read_length" in parameters["tool_options"]["cutadapt"]["hifi"] else "",
        min_adapter_length=lambda wildcards: " --overlap {0} ".format(parameters["tool_options"]["cutadapt"]["hifi"]["min_adapter_length"]) if "min_adapter_length" in parameters["tool_options"]["cutadapt"]["hifi"] else "",
        adapter_match_times=lambda wildcards: " --times {0}".format(parameters["tool_options"]["cutadapt"]["hifi"]["adapter_match_times"]) if "adapter_match_times" in parameters["tool_options"]["cutadapt"]["hifi"] else "",
        check_read_rc= lambda wildcards: " --rc " if ( ("check_read_rc" in parameters["tool_options"]["cutadapt"]["hifi"]) and  parameters["tool_options"]["cutadapt"]["hifi"]["check_read_rc"]) else "",
        discard_trimmed= lambda wildcards: " --discard-trimmed " if ( ("discard_trimmed" in parameters["tool_options"]["cutadapt"]["hifi"]) and parameters["tool_options"]["cutadapt"]["hifi"]["discard_trimmed"]) else "",
        forward_anywhere_adapters= lambda wildcards: (" -b " + " -b ".join(parameters["tool_options"]["cutadapt"]["hifi"]["forward_anywhere_adapter_list"])) if ("forward_anywhere_adapter_list" in parameters["tool_options"]["cutadapt"]["hifi"]) and parameters["tool_options"]["cutadapt"]["hifi"]["forward_anywhere_adapter_list"] else "",
        reverse_anywhere_adapters= lambda wildcards: (" -B " + " -B ".join(parameters["tool_options"]["cutadapt"]["hifi"]["reverse_anywhere_adapter_list"])) if ("reverse_anywhere_adapter_list" in parameters["tool_options"]["cutadapt"]["hifi"]) and parameters["tool_options"]["cutadapt"]["hifi"]["reverse_anywhere_adapter_list"] else "",
        forward_three_prime_adapters= lambda wildcards: (" -a " + " -a ".join(parameters["tool_options"]["cutadapt"]["hifi"]["forward_three_prime_adapter_list"])) if ("forward_three_prime_adapter_list" in parameters["tool_options"]["cutadapt"]["hifi"]) and parameters["tool_options"]["cutadapt"]["hifi"]["forward_three_prime_adapter_list"] else "",
        reverse_three_prime_adapters= lambda wildcards: (" -A " + " -A ".join(parameters["tool_options"]["cutadapt"]["hifi"]["reverse_three_prime_adapter_list"])) if ("reverse_three_prime_adapter_list" in parameters["tool_options"]["cutadapt"]["hifi"]) and parameters["tool_options"]["cutadapt"]["hifi"]["reverse_three_prime_adapter_list"] else "",
        forward_five_prime_adapters= lambda wildcards: (" -g " + " -g ".join(parameters["tool_options"]["cutadapt"]["hifi"]["forward_five_prime_adapter_list"])) if ("forward_five_prime_adapter_list" in parameters["tool_options"]["cutadapt"]["hifi"]) and parameters["tool_options"]["cutadapt"]["hifi"]["forward_five_prime_adapter_list"] else "",
        reverse_five_prime_adapters= lambda wildcards: (" -G " + " -G ".join(parameters["tool_options"]["cutadapt"]["hifi"]["reverse_five_prime_adapter_list"])) if ("reverse_five_prime_adapter_list" in parameters["tool_options"]["cutadapt"]["hifi"]) and parameters["tool_options"]["cutadapt"]["hifi"]["reverse_five_prime_adapter_list" ] else "",
    log:
        std=output_dict["log"] / "cutadapt_pacbio.hifi.{fileprefix}.log",
        cluster_log=output_dict["cluster_log"] / "cutadapt_pacbio.hifi.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "cutadapt_pacbio.hifi.{fileprefix}.cluster.log"
    benchmark:
        output_dict["benchmark"] /  "cutadapt_pacbio.hifi.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("cutadapt"),
        cpus=parameters["threads"]["cutadapt"],
        time=parameters["time"]["cutadapt"],
        mem=parameters["memory_mb"]["cutadapt"],
    threads:
        parameters["threads"]["cutadapt"]
    shell:
         " cutadapt -j {threads} {params.min_read_length} {params.error_rate} {params.min_adapter_length} "
         " {params.adapter_match_times} "
         " {params.forward_anywhere_adapters} {params.reverse_anywhere_adapters} "
         " {params.forward_five_prime_adapters} {params.reverse_five_prime_adapters} "
         " {params.forward_three_prime_adapters} {params.reverse_three_prime_adapters} "
         " {params.check_read_rc} {params.discard_trimmed} "
         "  -o {output.fastq} {input.fastq} > {output.stats} 2>{log.std}; "

#if "illumina" in data_types:
rule cutadapt_illumina:
    input:
        forward_fastq=lambda wildcards: output_dict["data"] / ("fastq/{0}/raw/{1}{2}{3}".format(wildcards.datatype,
                                                                                                wildcards.pairprefix,
                                                                                                input_forward_suffix_dict[wildcards.datatype],
                                                                                                config["fastq_extension"])),
        reverse_fastq=lambda wildcards: output_dict["data"] / ("fastq/{0}/raw/{1}{2}{3}".format(wildcards.datatype,
                                                                                                wildcards.pairprefix,
                                                                                                input_reverse_suffix_dict[wildcards.datatype],
                                                                                                config["fastq_extension"])),
    output:
        forward_fastq=output_dict["data"] / ("fastq/{datatype, hic|illumina}/filtered/{pairprefix}_1%s" % config["fastq_extension"]),
        reverse_fastq=output_dict["data"] / ("fastq/{datatype, hic|illumina}/filtered/{pairprefix}_2%s" % config["fastq_extension"]),
        stats=output_dict["data"] / "fastq/{datatype, hic|illumina}/filtered/{pairprefix}.cutadapt.stats"
    params:
        error_rate=lambda wildcards: "-e {0} ".format(parameters["tool_options"]["cutadapt"][wildcards.datatype]["error_rate"]) if "error_rate" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        min_read_length=lambda wildcards: " -m {0} ".format(parameters["tool_options"]["cutadapt"][wildcards.datatype]["min_read_length"]) if "min_read_length" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        min_adapter_length=lambda wildcards: " --overlap {0} ".format(parameters["tool_options"]["cutadapt"][wildcards.datatype]["min_adapter_length"]) if "min_adapter_length" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        adapter_match_times=lambda wildcards: " --times {0}".format(parameters["tool_options"]["cutadapt"][wildcards.datatype]["adapter_match_times"]) if "adapter_match_times" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        check_read_rc= lambda wildcards: " --rc " if ( ("check_read_rc" in parameters["tool_options"]["cutadapt"][wildcards.datatype]) and  parameters["tool_options"]["cutadapt"][wildcards.datatype]["check_read_rc"]) else "",
        discard_trimmed= lambda wildcards: " --discard-trimmed " if ( ("discard_trimmed" in parameters["tool_options"]["cutadapt"][wildcards.datatype]) and parameters["tool_options"]["cutadapt"][wildcards.datatype]["discard_trimmed"]) else "",
        forward_anywhere_adapters= lambda wildcards: (" -b " + " -b ".join(parameters["tool_options"]["cutadapt"][wildcards.datatype]["forward_anywhere_adapter_list"])) if "forward_anywhere_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        reverse_anywhere_adapters= lambda wildcards: (" -B " + " -B ".join(parameters["tool_options"]["cutadapt"][wildcards.datatype]["reverse_anywhere_adapter_list"])) if "reverse_anywhere_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        forward_three_prime_adapters= lambda wildcards: (" -a " + " -a ".join(parameters["tool_options"]["cutadapt"][wildcards.datatype]["forward_three_prime_adapter_list"])) if "forward_three_prime_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        reverse_three_prime_adapters= lambda wildcards: (" -A " + " -A ".join(parameters["tool_options"]["cutadapt"][wildcards.datatype]["reverse_three_prime_adapter_list"])) if "reverse_three_prime_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        forward_five_prime_adapters= lambda wildcards: (" -g " + " -g ".join(parameters["tool_options"]["cutadapt"][wildcards.datatype]["forward_five_prime_adapter_list"])) if "forward_five_prime_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
        reverse_five_prime_adapters= lambda wildcards: (" -G " + " -G ".join(parameters["tool_options"]["cutadapt"][wildcards.datatype]["reverse_five_prime_adapter_list"])) if "reverse_five_prime_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.datatype] else "",
    log:
        std=output_dict["log"] / "cutadapt_pacbio.{datatype}.{pairprefix}.log",
        #stats=log_dir_path / "{library_id}/no_cut.cutadapt.stats.log",
        cluster_log=output_dict["cluster_log"] / "cutadapt_pacbio.{datatype}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "cutadapt_pacbio.{datatype}.{pairprefix}.cluster.log"
    benchmark:
        output_dict["benchmark"] /  "cutadapt_pacbio.{datatype}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("cutadapt"),
        cpus=parameters["threads"]["cutadapt"],
        time=parameters["time"]["cutadapt"],
        mem=parameters["memory_mb"]["cutadapt"],
    threads:
        parameters["threads"]["cutadapt"]
    shell:
         " cutadapt -j {threads} {params.min_read_length} {params.error_rate} {params.min_adapter_length} "
         " {params.adapter_match_times} "
         " {params.forward_anywhere_adapters} {params.reverse_anywhere_adapters} "
         " {params.forward_five_prime_adapters} {params.reverse_five_prime_adapters} "
         " {params.forward_three_prime_adapters} {params.reverse_three_prime_adapters} "
         " {params.check_read_rc} {params.discard_trimmed} "
         " -o {output.forward_fastq} -p {output.reverse_fastq} "
         " {input.forward_fastq} {input.reverse_fastq} > {output.stats} 2>{log.std}; "
