
rule combine_se_reads:
    input:
        se_reads=lambda wildcards: expand(output_dict["data"] / ("fastq/%s/filtered/{fileprefix}%s" % (wildcards.datatype,
                                                                                                         config["fastq_extension"])),
                          fileprefix=input_file_prefix_dict[wildcards.datatype])
    output:
        combined_se_reads=output_dict["data"] / ("fastq/{datatype, hifi|nanopore|simplex|duplex}/combined/{datatype}.combined%s" % config["fastq_extension"])
    log:
        log=output_dict["log"] / "combine_long_reads.{datatype}.log",
        cluster_log=output_dict["cluster_log"] / "combine_long_reads.{datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "combine_long_reads.{datatype}.err"
    benchmark:
        output_dict["benchmark"] / "combine_long_reads.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" %config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("combine_long_reads"),
        cpus=parameters["threads"]["combine_long_reads"],
        time=parameters["time"]["combine_long_reads"],
        mem=parameters["memory_mb"]["combine_long_reads"],
    threads: parameters["threads"]["combine_long_reads"]
    shell:
        "cat {input.se_reads} > {output.combined_se_reads} 2>{log.log}; "

rule combine_paired_reads:
    input:
        pe_forward_reads=lambda wildcards: expand(output_dict["data"] / ("fastq/%s/filtered/{pairprefix}_1%s" % (wildcards.datatype,
                                                                                                                 config["fastq_extension"])),
                                                  pairprefix=input_pairprefix_dict[wildcards.datatype]),
        pe_reverse_reads=lambda wildcards: expand(output_dict["data"] / ("fastq/%s/filtered/{pairprefix}_2%s" % (wildcards.datatype,
                                                                                                                 config["fastq_extension"])),
                                                  pairprefix=input_pairprefix_dict[wildcards.datatype]),
    output:
        combined_pe_forward_reads=output_dict["data"] / ("fastq/{datatype, hic|illumina}/combined/{datatype}.combined_1%s" % config["fastq_extension"]),
        combined_pe_reverse_reads=output_dict["data"] / ("fastq/{datatype, hic|illumina}/combined/{datatype}.combined_2%s" % config["fastq_extension"]),
    log:
        log=output_dict["log"] / "combine_long_reads.{datatype}.log",
        cluster_log=output_dict["cluster_log"] / "combine_long_reads.{datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "combine_long_reads.{datatype}.err"
    benchmark:
        output_dict["benchmark"] / "combine_long_reads.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" %config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("combine_long_reads"),
        cpus=parameters["threads"]["combine_long_reads"],
        time=parameters["time"]["combine_long_reads"],
        mem=parameters["memory_mb"]["combine_long_reads"],
    threads: parameters["threads"]["combine_long_reads"]
    shell:
        "cat {input.pe_forward_reads} > {output.combined_pe_forward_reads} 2>{log.log}; "
        "cat {input.pe_reverse_reads} > {output.combined_pe_reverse_reads} 2>>{log.log}; "

