localrules: merge_tadbit_stats
ruleorder: merge_tadbit_stats > tadbit

rule tadbit:
    input:
        forward_fastq=config["out_dir"] / ("data/{datatype}/{stage}/{pairprefix}%s%s" % (config["fwd_fastq_sfx"], config["fastq_ext"])),
        reverse_fastq=config["out_dir"] / ("data/{datatype}/{stage}/{pairprefix}%s%s" % (config["rev_fastq_sfx"], config["fastq_ext"])),
    output:
        stats=config["out_dir"] / "qc/tadbit/{datatype, hic}/{stage}/{pairprefix}.data.stats" ,
        #stats=merged_raw_fastqc_dir_path / "{library_id}/{library_id}.raw.fast{}qc.stats"
    params:
        enzyme_list=",".join(config["hic_enzyme_dict"][config["hic_enzyme_set"]] if config["custom_enzyme_set"] is None else config["custom_enzyme_set"]),
        read_number=parameters["tool_options"]["tadbit"]["hic"]["read_number"]
    log:
        forward_tadbit=config["out_dir"] / "log/tadbit.{stage}.{datatype}.{pairprefix}.forward.log",
        reverse_tadbit=config["out_dir"] / "log/tadbit.{stage}.{datatype}.{pairprefix}.reverse.log",
        combine=config["out_dir"] / "log/tadbit.{stage}.{datatype}.{pairprefix}.combine.log",
        cluster_log=config["out_dir"] / "log/tadbit.{stage}.{datatype}.{pairprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/tadbit.{stage}.{datatype}.{pairprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/tadbit.{stage}.{datatype}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["tadbit"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["tadbit"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
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
        "     -e {params.enzyme_list} -p ${{OUTPUT_PREFIX}}.forward > {log.forward_tadbit} 2>&1; "
        " workflow/scripts/hic_qc/count_ligation_site_metrics.py  -s 'reverse_' -n {params.read_number} -f {input.reverse_fastq} "
        "     -e {params.enzyme_list} -p ${{OUTPUT_PREFIX}}.reverse > {log.reverse_tadbit} 2>&1; "
        " paste ${{OUTPUT_PREFIX}}.forward.stats <(cut -f 3 ${{OUTPUT_PREFIX}}.reverse.stats) > ${{OUTPUT_PREFIX}}.stats 2>{log.combine}; "

rule merge_tadbit_stats:
    input:
        stats=lambda wildcards: expand(config["out_dir"] / "qc/tadbit/{datatype}/{stage}/{pairprefix}.data.stats",
                                       datatype=[wildcards.datatype,],
                                       pairprefix=config["data"][wildcards.datatype]["pair_prefix_list"],
                                       allow_missing=True)
    output:
        stats=config["out_dir"] / "qc/tadbit/{datatype, hic}/{stage}/{genome_prefix}.tadbit.stats" ,
    log:
        head=config["out_dir"] / "log/merge_tadbit_stats.{stage}.{datatype}.{genome_prefix}.head.log",
        tail=config["out_dir"] / "log/merge_tadbit_stats.{stage}.{datatype}.{genome_prefix}.tail.log",
        cluster_log=config["out_dir"] / "log/tadbit.{stage}.{datatype}.{genome_prefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/tadbit.{stage}.{datatype}.{genome_prefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/tadbit.{stage}.{datatype}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("merge_tadbit_stats"),
        cpus=parameters["threads"]["merge_tadbit_stats"],
        time=parameters["time"]["merge_tadbit_stats"],
        mem=parameters["memory_mb"]["merge_tadbit_stats"],
    threads:
        parameters["threads"]["merge_tadbit_stats"]
    shell:
        " > {log.tail}; "
        " INPUT_FILES_ARR=({input.stats}); "
        " head -n 1 ${{INPUT_FILES_ARR[0]}} | sed 's/^#/#genome_prefix\tpair_prefix\t/' > {output.stats} 2>{log.head}; "
        " for STAT_FILE in {input.stats};"
        " do "
        "     PAIR_PREFIX=`basename ${{STAT_FILE}}`; "
        "     PAIR_PREFIX=${{PAIR_PREFIX%.stats}}; "
        "     tail -n +2 ${{STAT_FILE}} | "
        "     awk -v PAIR_PREFIX=${{PAIR_PREFIX}} '{{print \"{wildcards.genome_prefix}\t\"PAIR_PREFIX\"\t\"$0}}' >> {output.stats} 2>>{log.tail}; "
        " done "