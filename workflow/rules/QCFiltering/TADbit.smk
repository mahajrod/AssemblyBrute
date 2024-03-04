localrules: merge_tadbit_stats

rule tadbit:
    input:
        #fastq_dir=rules.create_fastq_links.output,
        forward_fastq=output_dict["data"] / ("fastq/hic/raw/{pairprefix}%s%s" % (input_forward_suffix_dict["hic"], 
                                                                                     config["fastq_extension"])),
        reverse_fastq=output_dict["data"] / ("fastq/hic/raw/{pairprefix}%s%s" % (input_reverse_suffix_dict["hic"], 
                                                                                     config["fastq_extension"])),
    output:
        stats=output_dict["qc"] / "tadbit/hic/raw/{pairprefix}.stats" ,
        #stats=merged_raw_fastqc_dir_path / "{library_id}/{library_id}.raw.fast{}qc.stats"
    params:
        enzyme_list=",".join(config["hic_enzyme_dict"][config["hic_enzyme_set"]] if config["custom_enzyme_set"] is None else config["custom_enzyme_set"]),
        read_number=parameters["tool_options"]["tadbit"]["hic"]["read_number"]
    log:
        forward_tadbit=output_dict["log"]/ "tadbit.raw.{pairprefix}.forward.log",
        reverse_tadbit=output_dict["log"]/ "tadbit.raw.{pairprefix}.reverse.log",
        combine=output_dict["log"]/ "tadbit.raw.{pairprefix}.combine.log",
        #stats=log_dir_path / "{library_id}/fastqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"] / "tadbit.raw.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "tadbit.raw.{pairprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "tadbit.raw.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["tadbit"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["tadbit"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("tadbit"),
        cpus=parameters["threads"]["tadbit"],
        time=parameters["time"]["tadbit"],
        mem=parameters["memory_mb"]["tadbit"],
    threads:
        parameters["threads"]["tadbit"]
    shell:
        " OUTPUT_PREFIX={output.stats}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.stats}}; "
        " workflow/scripts/hic_qc/count_ligation_site_metrics.py  -s 'forward_' -n {params.read_number} -f {input.forward_fastq} "
        " -e {params.enzyme_list} -p ${{OUTPUT_PREFIX}}.forward > {log.forward_tadbit} 2>&1; "
        " workflow/scripts/hic_qc/count_ligation_site_metrics.py  -s 'reverse_' -n {params.read_number} -f {input.reverse_fastq} "
        " -e {params.enzyme_list} -p ${{OUTPUT_PREFIX}}.reverse > {log.reverse_tadbit} 2>&1; "
        " paste ${{OUTPUT_PREFIX}}.forward.stats <(cut -f 3 ${{OUTPUT_PREFIX}}.reverse.stats) > ${{OUTPUT_PREFIX}} 2>{log.combine}; "

rule merge_tadbit_stats:
    input:
        stats=expand(output_dict["qc"] / "tadbit/hic/raw/{pairprefix}.stats", pairprefix=input_pairprefix_dict["hic"])
    output:
        stats=output_dict["qc"] / "tadbit/hic/raw/{genome_prefix}.tadbit.stats" ,
        #stats=merged_raw_fastqc_dir_path / "{library_id}/{library_id}.raw.fast{}qc.stats"
    params:
        header_file=expand(output_dict["qc"] / "tadbit/hic/raw/{pairprefix}.stats", pairprefix=input_pairprefix_dict["hic"])[0]
    log:
        head=output_dict["log"]/ "merge_tadbit_stats.raw.{genome_prefix}.head.log",
        tail=output_dict["log"]/ "merge_tadbit_stats.raw.{genome_prefix}.tail.log",
        #stats=log_dir_path / "{library_id}/fastqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"] / "tadbit.raw.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "tadbit.raw.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "tadbit.raw.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("merge_tadbit_stats"),
        cpus=parameters["threads"]["merge_tadbit_stats"],
        time=parameters["time"]["merge_tadbit_stats"],
        mem=parameters["memory_mb"]["merge_tadbit_stats"],
    threads:
        parameters["threads"]["merge_tadbit_stats"]
    shell:
        " > {log.tail}; "
        " head -n 1 {params.header_file} > {output.stats} 2>{log.head}; "
        " for STAT_FILE in {input.stats};"
        "   do "
        "   tail -n +2 ${{STAT_FILE}} >> {output.stats} 2>>{log.tail}; "
        "   done "