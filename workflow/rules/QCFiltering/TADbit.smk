localrules: merge_tadbit_stats
ruleorder: merge_tadbit_stats > tadbit

rule tadbit:
    input:
        #fastq_dir=rules.create_fastq_links.output,
        forward_fastq=lambda wildcards: output_dict["data"] / ("fastq/{0}/raw/{1}{2}{3}".format(wildcards.datatype,
                                                                                              wildcards.pairprefix,
                                                                                              input_forward_suffix_dict[wildcards.datatype],
                                                                                              config["fastq_extension"])),
        reverse_fastq=lambda wildcards: output_dict["data"] / ("fastq/{0}/raw/{1}{2}{3}".format(wildcards.datatype,
                                                                                              wildcards.pairprefix,
                                                                                              input_reverse_suffix_dict[wildcards.datatype],
                                                                                              config["fastq_extension"])),
    output:
        stats=output_dict["qc"] / "tadbit/{datatype, hic}/raw/{pairprefix}.stats" ,
        #stats=merged_raw_fastqc_dir_path / "{library_id}/{library_id}.raw.fast{}qc.stats"
    params:
        enzyme_list=",".join(config["hic_enzyme_dict"][config["hic_enzyme_set"]] if config["custom_enzyme_set"] is None else config["custom_enzyme_set"]),
        read_number=parameters["tool_options"]["tadbit"]["hic"]["read_number"]
    log:
        forward_tadbit=output_dict["log"]/ "tadbit.raw.{datatype}.{pairprefix}.forward.log",
        reverse_tadbit=output_dict["log"]/ "tadbit.raw.{datatype}.{pairprefix}.reverse.log",
        combine=output_dict["log"]/ "tadbit.raw.{datatype}.{pairprefix}.combine.log",
        #stats=log_dir_path / "{library_id}/fastqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"] / "tadbit.raw.{datatype}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "tadbit.raw.{datatype}.{pairprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "tadbit.raw.{datatype}.{pairprefix}.benchmark.txt"
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
        " paste ${{OUTPUT_PREFIX}}.forward.stats <(cut -f 3 ${{OUTPUT_PREFIX}}.reverse.stats) > ${{OUTPUT_PREFIX}}.stats 2>{log.combine}; "

rule merge_tadbit_stats:
    input:
        stats=lambda wildcards: expand(output_dict["qc"] / "tadbit/%s/raw/{pairprefix}.stats" % wildcards.datatype,
                                       pairprefix=input_pairprefix_dict[wildcards.datatype])
    output:
        stats=output_dict["qc"] / "tadbit/{datatype, hic}/raw/{genome_prefix}.tadbit.stats" ,
        #stats=merged_raw_fastqc_dir_path / "{library_id}/{library_id}.raw.fast{}qc.stats"
    #params:
    #    header_file=expand(output_dict["qc"] / "tadbit/{datatype}/raw/{pairprefix}.stats", pairprefix=input_pairprefix_dict["hic"])[0]
    log:
        head=output_dict["log"]/ "merge_tadbit_stats.raw.{datatype}.{genome_prefix}.head.log",
        tail=output_dict["log"]/ "merge_tadbit_stats.raw.{datatype}.{genome_prefix}.tail.log",
        #stats=log_dir_path / "{library_id}/fastqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"] / "tadbit.raw.{datatype}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "tadbit.raw.{datatype}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "tadbit.raw.{datatype}.{genome_prefix}.benchmark.txt"
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
        " INPUT_FILES_ARR=({input.stats}); "
        " head -n 1 ${{ARR[0]}} | sed 's/^#/#genome_prefix\tpair_prefix\t/' > {output.stats} 2>{log.head}; "
        " for STAT_FILE in {input.stats};"
        "   do "
        "   PAIR_PREFIX=`basename ${{STAT_FILE}}`; "
        "   PAIR_PREFIX=${{PAIR_PREFIX%.stats}}; "
        "   tail -n +2 ${{STAT_FILE}} | "
        "   awk -v PAIR_PREFIX=${{PAIR_PREFIX}} '{{print \"{wildcards.genome_prefix}\t\"PAIR_PREFIX\"\t\"$0}}' >> {output.stats} 2>>{log.tail}; "
        "   done "