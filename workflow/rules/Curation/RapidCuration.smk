localrules: create_curation_input_links

rule create_curation_input_links: #
    input:
        fasta=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta" % (stage_dict["curation"]["prev_stage"],
                                                                                                   stage_dict["curation"]["prev_stage"])),
        fai=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta.fai" % (stage_dict["curation"]["prev_stage"],
                                                                                                       stage_dict["curation"]["prev_stage"]))
    output:
        fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.fasta",
        fai=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.fasta.fai"
    log:
        ln1=output_dict["log"]  / "create_curation_input_links.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.ln1.log",
        ln2=output_dict["log"]  / "create_curation_input_links.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.ln2.log",
        cluster_log=output_dict["cluster_log"] / "create_curation_input_links.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_curation_input_links.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_curation_input_links.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["create_curation_input_links"],
        time=parameters["time"]["create_curation_input_links"],
        mem=parameters["memory_mb"]["create_curation_input_links"]
    threads: parameters["threads"]["create_curation_input_links"]

    shell:
        " ln -sf `realpath -s {input.fasta}` {output.fasta} > {log.ln1} 2>&1; "
        " ln -sf `realpath -s {input.fai}` {output.fai} > {log.ln2} 2>&1; "
