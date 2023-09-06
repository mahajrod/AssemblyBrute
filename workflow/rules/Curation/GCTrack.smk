
ruleorder: create_gc_track > create_bedgraph_track
rule create_gc_track: #
    input:
        fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.fasta"
    output:
        gap_bedgraph=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/{seq_type}/{genome_prefix}.input.{haplotype}.gc.win{window}.step{step}.track.bedgraph",
    log:
        gc=output_dict["log"]  / "create_gc_track.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.win{window}.step{step}.gc.log",
        cluster_log=output_dict["cluster_log"] / "create_gc_track.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.win{window}.step{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_gc_track.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.win{window}.step{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_gc_track.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.win{window}.step{step}.benchmark.txt"
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
        " -w {wildcards.window} -s {wildcards.step} -o {output.gap_bedgraph} 2>{log.gc}; "
