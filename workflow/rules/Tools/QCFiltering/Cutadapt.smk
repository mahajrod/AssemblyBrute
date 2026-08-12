
rule cutadapt:
    input:
        fastq=config["out_dir"] / ("data/{pacbio_datatype}/raw/{fileprefix}%s" % config["fastq_ext"])
    output:
        fastq=config["out_dir"] / ("data/{pacbio_datatype}/filtered/{fileprefix}%s" % config["fastq_ext"]),
        stats=config["out_dir"] / "data/{pacbio_datatype}/filtered/{fileprefix}.cutadapt.stats"
    params:
        error_rate=lambda wildcards: "-e {0} ".format(parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["error_rate"]) if "error_rate" in parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype] else "",
        min_read_length=lambda wildcards: " -m {0} ".format(parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["min_read_length"]) if "min_read_length" in parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype] else "",
        min_adapter_length=lambda wildcards: " --overlap {0} ".format(parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["min_adapter_length"]) if "min_adapter_length" in parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype] else "",
        adapter_match_times=lambda wildcards: " --times {0}".format(parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["adapter_match_times"]) if "adapter_match_times" in parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype] else "",
        check_read_rc= lambda wildcards: " --rc " if ( ("check_read_rc" in parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]) and  parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["check_read_rc"]) else "",
        discard_trimmed= lambda wildcards: " --discard-trimmed " if ( ("discard_trimmed" in parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]) and parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["discard_trimmed"]) else "",
        forward_anywhere_adapters= lambda wildcards: (" -b " + " -b ".join(parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["forward_anywhere_adapter_list"])) if ("forward_anywhere_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]) and parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["forward_anywhere_adapter_list"] else "",
        reverse_anywhere_adapters= lambda wildcards: (" -B " + " -B ".join(parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["reverse_anywhere_adapter_list"])) if ("reverse_anywhere_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]) and parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["reverse_anywhere_adapter_list"] else "",
        forward_three_prime_adapters= lambda wildcards: (" -a " + " -a ".join(parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["forward_three_prime_adapter_list"])) if ("forward_three_prime_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]) and parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["forward_three_prime_adapter_list"] else "",
        reverse_three_prime_adapters= lambda wildcards: (" -A " + " -A ".join(parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["reverse_three_prime_adapter_list"])) if ("reverse_three_prime_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]) and parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["reverse_three_prime_adapter_list"] else "",
        forward_five_prime_adapters= lambda wildcards: (" -g " + " -g ".join(parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["forward_five_prime_adapter_list"])) if ("forward_five_prime_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]) and parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["forward_five_prime_adapter_list"] else "",
        reverse_five_prime_adapters= lambda wildcards: (" -G " + " -G ".join(parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["reverse_five_prime_adapter_list"])) if ("reverse_five_prime_adapter_list" in parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]) and parameters["tool_options"]["cutadapt"][wildcards.pacbio_datatype]["reverse_five_prime_adapter_list" ] else "",
    log:
        std=config["out_dir"] / "log/cutadapt.{pacbio_datatype}.{fileprefix}.log",
        cluster_log=config["out_dir"] / "log/cutadapt.hifi.{pacbio_datatype}.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/cutadapt.hifi.{pacbio_datatype}.{fileprefix}.cluster.log"
    benchmark:
        config["out_dir"] / "log/cutadapt.hifi.{pacbio_datatype}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("cutadapt"),
        cpus=parameters["threads"]["cutadapt"],
        time=parameters["time"]["cutadapt"],
        mem=parameters["memory_mb"]["cutadapt"],
    threads:
        parameters["threads"]["cutadapt"]
    shell:
         " cutadapt -j {threads} {params.min_read_length} {params.error_rate} {params.min_adapter_length} "
         "     {params.adapter_match_times} {params.forward_anywhere_adapters} {params.reverse_anywhere_adapters} "
         "     {params.forward_five_prime_adapters} {params.reverse_five_prime_adapters} "
         "     {params.forward_three_prime_adapters} {params.reverse_three_prime_adapters} "
         "     {params.check_read_rc} {params.discard_trimmed} -o {output.fastq} {input.fastq} > {output.stats} 2>{log.std}; "

use rule cutadapt as cutadapt_track_data with:
    input:
        fastq=config["out_dir"] / ("ext_track_data/{pacbio_datatype}/{track_name}/raw/{fileprefix}%s" % config["fastq_ext"])
    output:
        fastq=config["out_dir"] / ("ext_track_data/{pacbio_datatype}/{track_name}/filtered/{fileprefix}%s" % config["fastq_ext"]),
        stats=config["out_dir"] / "ext_track_data/{pacbio_datatype}/{track_name}/filtered/{fileprefix}.cutadapt.stats"
    log:
        std=config["out_dir"] / "log/cutadapt_track_data.{pacbio_datatype}.{track_name}.{fileprefix}.log",
        cluster_log=config["out_dir"] / "log/cutadapt_track_data.hifi.{pacbio_datatype}.{track_name}.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/cutadapt_track_data.hifi.{pacbio_datatype}.{track_name}.{fileprefix}.cluster.log"
    benchmark:
        config["out_dir"] / "log/cutadapt_track_data.hifi.{pacbio_datatype}.{track_name}.{fileprefix}.benchmark.txt"