
rule combine_se_reads:
    input:
        se_reads=lambda wildcards: expand(config["out_dir"] / ("data/%s/final/{fileprefix}%s" % (wildcards.se_datatype,
                                                                                                       config["data"][wildcards.se_datatype]["conv_ext"])),
                          fileprefix=config["data"][wildcards.se_datatype]["conv_file_prefix_list"])
    output:
        combined_se_reads=config["out_dir"] / ("data/{se_datatype}/combined/{se_datatype}.combined%s" % config["fastq_ext"])
    params:
        number_of_files=lambda wildcards: len(config["data"][wildcards.se_datatype]["conv_file_prefix_list"])
    log:
        log=config["out_dir"] / "log/combine_long_reads.{se_datatype}.log",
        cluster_log=config["out_dir"] / "log/combine_long_reads.{se_datatype}.cluster.log",
        cluster_err=config["out_dir"] / "log/combine_long_reads.{se_datatype}.err"

    benchmark:
        config["out_dir"] / "log/combine_long_reads.{se_datatype}.benchmark.txt"
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
        " if [[ {params.number_of_files} == 1 ]]; "
        " then "
        "     ln -sf ../final/`basename {input.se_reads}` {output.combined_se_reads} > {log.log} 2>&1; "
        " else "
        "     cat {input.se_reads} > {output.combined_se_reads} 2>{log.log}; "
        " fi "

rule combine_paired_reads:
    input:
        pe_forward_reads=lambda wildcards: expand(config["out_dir"] / ("data/%s/filtered/{pairprefix}%s%s" % (wildcards.pe_datatype,
                                                                                                                    config["fwd_fastq_sfx"],
                                                                                                                    config["fastq_ext"])),
                                                  pairprefix=config["data"][wildcards.pe_datatype]["conv_pair_prefix_list"]),
        pe_reverse_reads=lambda wildcards: expand(config["out_dir"] / ("data/%s/filtered/{pairprefix}%s%s" % (wildcards.pe_datatype,
                                                                                                                        config["fwd_fastq_sfx"],
                                                                                                                        config["fastq_ext"])),
                                                  pairprefix=config["data"][wildcards.pe_datatype]["conv_pair_prefix_list"]),
    output:
        combined_pe_forward_reads=config["out_dir"] / ("data/{pe_datatype}/combined/{pe_datatype}.combined%s%s" % (config["fwd_fastq_sfx"], config["fastq_ext"])),
        combined_pe_reverse_reads=config["out_dir"] / ("data/{pe_datatype}/combined/{pe_datatype}.combined%s%s" % (config["rev_fastq_sfx"], config["fastq_ext"])),
    log:
        log=config["out_dir"] / "log/combine_long_reads.{pe_datatype}.log",
        cluster_log=config["out_dir"] / "log/combine_long_reads.{pe_datatype}.cluster.log",
        cluster_err=config["out_dir"] / "log/combine_long_reads.{pe_datatype}.err"
    benchmark:
        config["out_dir"] / "log/combine_long_reads.{pe_datatype}.benchmark.txt"
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
