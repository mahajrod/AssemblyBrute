

rule create_gap_track: #
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta"
    output:
        gap_bed=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/{track_type, gap}/{haplotype, [^.]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.{track_type}.track.bed",
    log:
        seqtk=output_dict["log"]  / "create_gap_track.{assembly_stage}..{parameters}.{track_type}.{genome_prefix}.{haplotype}.seqtk.log",
        cluster_log=output_dict["cluster_log"] / "create_gap_track.{assembly_stage}..{parameters}.{track_type}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_gap_track.{assembly_stage}..{parameters}.{track_type}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_gap_track.{assembly_stage}..{parameters}.{track_type}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("create_gap_track"),
        cpus=parameters["threads"]["create_gap_track"],
        time=parameters["time"]["create_gap_track"],
        mem=parameters["memory_mb"]["create_gap_track"]
    threads: parameters["threads"]["create_gap_track"]

    shell:
        " seqtk cutN -n 1 -g  {input.fasta} > {output.gap_bed} 2>{log.seqtk} "
