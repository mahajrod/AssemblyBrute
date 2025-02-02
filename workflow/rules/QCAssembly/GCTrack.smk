
ruleorder: create_gc_track > create_bedgraph_track
rule create_gc_track: #
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta"
    output:
        gc_bedgraph=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/tracks/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^.]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.{track_type, gc}.win{window, [0-9]+}.step{step, [0-9]+}.track.bedgraph",
    log:
        gc=output_dict["log"]  / "create_gc_track.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.win{window}.step{step}.gc.log",
        cluster_log=output_dict["cluster_log"] / "create_gc_track.{assembly_stage}.{track_type}.{parameters}.{genome_prefix}.{haplotype}.win{window}.step{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_gc_track.{assembly_stage}.{track_type}.{parameters}.{genome_prefix}.{haplotype}.win{window}.step{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_gc_track.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.win{window}.step{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("create_gc_track"),
        cpus=parameters["threads"]["create_gc_track"],
        time=parameters["time"]["create_gc_track"],
        mem=parameters["memory_mb"]["create_gc_track"]
    threads: parameters["threads"]["create_gc_track"]

    shell:
        " workflow/scripts/curation/count_gc_in_windows.py -i {input.fasta} "
        " -w {wildcards.window} -s {wildcards.step} -o {output.gc_bedgraph} 2>{log.gc}; "
