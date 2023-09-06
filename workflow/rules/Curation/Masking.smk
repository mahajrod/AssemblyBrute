
rule maskfasta: #
    input:
        fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.fasta",
        trf_bed=rules.trf.output.bed,
        windowmasker_bed=rules.windowmasker.output.bed
    output:
        masked_fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/{seq_type}/{genome_prefix}.input.{haplotype}.softmasked.fasta",
        merged_bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/{seq_type}/{genome_prefix}.input.{haplotype}.repeats.track.bed",
    log:
        cat=output_dict["log"]  / "maskfasta.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cat.log",
        sort=output_dict["log"]  / "maskfasta.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.sort.log",
        merge=output_dict["log"]  / "maskfasta.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.merge.log",
        maskfasta=output_dict["log"]  / "maskfasta.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.maskfasta.log",
        cluster_log=output_dict["cluster_log"] / "maskfasta.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "maskfasta.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "maskfasta.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("maskfasta"),
        cpus=parameters["threads"]["maskfasta"] ,
        time=parameters["time"]["maskfasta"],
        mem=parameters["memory_mb"]["maskfasta"]
    threads: parameters["threads"]["maskfasta"]
    shell:
        " cat {input.trf_bed} {input.windowmasker_bed} 2>{log.cat} | "
        " sort -k1,1V -k2,2n -k3,3n 2>{log.sort} | "
        " bedtools merge -i stdin > {output.merged_bed} 2>{log.merge}; "
        " bedtools maskfasta -soft -fi {input.fasta} -bed {output.merged_bed} -fo {output.masked_fasta} > {log.maskfasta} 2>&1; "
