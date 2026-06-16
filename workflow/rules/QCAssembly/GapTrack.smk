

rule create_gap_track: #
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta"
    output:
        gap_bed=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/{track_type, gap}/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.track.bed",
        gap_bedgraph=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/{track_type, gap}/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.track.bedgraph",
        gap_bedgraph_alias=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type, gap}.track.bedgraph",

    log:
        seqtk=output_dict["log"]  / "create_gap_track.{assembly_stage}..{parameters}.{track_type}.{genome_prefix}.{haplotype}.seqtk.log",
        awk=output_dict["log"] / "create_gap_track.{assembly_stage}..{parameters}.{track_type}.{genome_prefix}.{haplotype}.awk.log",
        cp=output_dict["log"] / "create_gap_track.{assembly_stage}..{parameters}.{track_type}.{genome_prefix}.{haplotype}.cp.log",
        cluster_log=output_dict["cluster_log"] / "create_gap_track.{assembly_stage}..{parameters}.{track_type}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_gap_track.{assembly_stage}..{parameters}.{track_type}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_gap_track.{assembly_stage}..{parameters}.{track_type}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_gap_track"),
        cpus=parameters["threads"]["create_gap_track"],
        time=parameters["time"]["create_gap_track"],
        mem=parameters["memory_mb"]["create_gap_track"]
    threads: parameters["threads"]["create_gap_track"]

    shell:
        " seqtk cutN -n 1 -g  {input.fasta} > {output.gap_bed} 2>{log.seqtk}; "
        " awk '{{printf \"%s\\t%i\\n\",$0,1}}' {output.gap_bed} > {output.gap_bedgraph} 2>{log.awk}; "
        " cp {output.gap_bedgraph} {output.gap_bedgraph_alias} > {log.cp} 2>&1; "
