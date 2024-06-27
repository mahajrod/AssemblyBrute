
rule krater_from_histo:
    input:
        histo=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.{kmer_tool}.histo"
    output:

        summary=output_dict["kmer"] / "{datatype}/{stage}/krater/{datatype}.{stage}.{kmer_length}.{kmer_tool}/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.histo.stats",
        summary_alias=output_dict["kmer"] / "{datatype}/{stage}/krater/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.krater.parameters",
        local_maximums=output_dict["kmer"] / "{datatype}/{stage}/krater/{datatype}.{stage}.{kmer_length}.{kmer_tool}/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.local_maximums",
        local_minimums=output_dict["kmer"] / "{datatype}/{stage}/krater/{datatype}.{stage}.{kmer_length}.{kmer_tool}/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.local_minimums",
        logscale_png=output_dict["kmer"] / "{datatype}/{stage}/krater/{datatype}.{stage}.{kmer_length}.{kmer_tool}/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.logscale.png",
        normal_scale_png=output_dict["kmer"] / "{datatype}/{stage}/krater/{datatype}.{stage}.{kmer_length}.{kmer_tool}/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.no_logscale.png",
        both_png=output_dict["kmer"] / "{datatype}/{stage}/krater/{datatype}.{stage}.{kmer_length}.{kmer_tool}/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.png",
        both_peaks_and_gaps_png=output_dict["kmer"] / "{datatype}/{stage}/krater/{datatype}.{stage}.{kmer_length}.{kmer_tool}/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.peaks_and_gaps.png"

    params:
        #max_coverage=lambda wildcards: parameters["tool_options"][wildcards.kmer_tool][wildcards.datatype]["max_coverage"],
        low_limit=10, # TODO: add as option in config
        high_limit=150, # TODO: add as option in config
    log:
        std=output_dict["log"] / "krater_from_histo.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.log",
        cluster_log=output_dict["cluster_log"] / "krater_from_histo.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "krater_from_histo.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "krater_from_histo.{datatype}.{stage}.{kmer_length}.{kmer_tool}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("krater_from_histo"),
        cpus=parameters["threads"]["krater_from_histo"],
        time=parameters["time"]["krater_from_histo"],
        mem=parameters["memory_mb"]["krater_from_histo"],
    threads:
        parameters["threads"]["krater_from_histo"]
    shell:
         " OUTPUT_PREFIX={output.summary}; "
         " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.histo.stats}}; "
         " ~/Dropbox/KrATER/scripts/jf/draw_kmer_distribution_from_histo.py -i {input.histo} "
         " -a {wildcards.genome_prefix} -o ${{OUTPUT_PREFIX}} -w {params.low_limit} -g {params.high_limit} "
         " -m {wildcards.kmer_length} -d -n --dont_show_genome_size_on_plot > {log.std} 2>&1; "
         " cp -f {output.summary} {output.summary_alias} >> {log.std} 2>&1; " # -m {params.max_coverage}

