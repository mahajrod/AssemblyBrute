ruleorder: trimmomatic_pe> cutadapt > porechop_abi

rule trimmomatic_pe:
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
        forward_fastq=output_dict["data"] / ("fastq/{datatype, hic|illumina}/filtered/{pairprefix, [^/]+}_1%s" % config["fastq_extension"]),
        forward_se_fastq=output_dict["data"] / ("fastq/{datatype, hic|illumina}/filtered/{pairprefix, [^/]+}_1.se%s" % config["fastq_extension"]),
        reverse_fastq=output_dict["data"] / ("fastq/{datatype, hic|illumina}/filtered/{pairprefix, [^/]+}_2%s" % config["fastq_extension"]),
        reverse_se_fastq=output_dict["data"] / ("fastq/{datatype, hic|illumina}/filtered/{pairprefix, [^/]+}_2.se%s" % config["fastq_extension"]),
        stats=output_dict["data"] / "fastq/{datatype, hic|illumina}/filtered/{pairprefix, [^/]+}.trimmomatic.stats"
    params:
        min_read_length=lambda wildcards: parameters["tool_options"]["trimmomatic"][wildcards.datatype]["min_read_length"],
        sliding_window_size=lambda wildcards: parameters["tool_options"]["trimmomatic"][wildcards.datatype]["sliding_window_size"],
        sliding_window_quality=lambda wildcards: parameters["tool_options"]["trimmomatic"][wildcards.datatype]["sliding_window_quality"],
        illumina_clip=lambda wildcards: parameters["tool_options"]["trimmomatic"][wildcards.datatype]["illumina_clip"],
        adapter_file=lambda wildcards: parameters["tool_options"]["trimmomatic"][wildcards.datatype]["adapter_file"],
    log:
        std=output_dict["log"] / "trimmomatic_pe.{datatype}.{pairprefix}.log",
        #stats=log_dir_path / "{library_id}/no_cut.cutadapt.stats.log",
        cluster_log=output_dict["cluster_log"] / "trimmomatic_pe.{datatype}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "trimmomatic_pe.{datatype}.{pairprefix}.cluster.log"
    benchmark:
        output_dict["benchmark"] /  "trimmomatic_pe.{datatype}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("trimmomatic_pe"),
        cpus=parameters["threads"]["trimmomatic_pe"],
        time=parameters["time"]["trimmomatic_pe"],
        mem=parameters["memory_mb"]["trimmomatic_pe"],
    threads:
        parameters["threads"]["trimmomatic_pe"]
    shell:
         "trimmomatic PE -threads {threads} -phred33 {input.forward_fastq} {input.reverse_fastq} "
         " {output.forward_fastq} {output.forward_se_fastq} {output.reverse_fastq} {output.reverse_se_fastq} "
         " ILLUMINACLIP:{params.adapter_file}:{params.illumina_clip} "
         " SLIDINGWINDOW:{params.sliding_window_size}:{params.sliding_window_quality} "
         " MINLEN:{params.min_read_length} > {output.stats} 2>{log.std}; "
