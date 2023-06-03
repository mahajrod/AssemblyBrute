

rule create_gap_track: #
    input:
        fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.fasta"
    output:
        gap_bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.gap.bed",
        gap_bedgraph=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.gap.bedgraph"
    log:
        seqtk=output_dict["log"]  / "create_gap_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.seqtk.log",
        tee=output_dict["log"]  / "create_gap_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.tee.log",
        awk=output_dict["log"]  / "create_gap_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.awk.log",
        cluster_log=output_dict["cluster_log"] / "create_gap_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_gap_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_gap_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["create_gap_track"],
        time=parameters["time"]["create_gap_track"],
        mem=parameters["memory_mb"]["create_gap_track"]
    threads: parameters["threads"]["create_gap_track"]

    shell:
        " seqtk cutN -n 1 -g  {input.fasta} 2>{log.seqtk} | tee {output.gap_bed} 2>{log.tee} | "
        " awk '{{print $0\"\t\"($3-$2)}}' > {output.gap_bedgraph} 2>{log.awk}; "
