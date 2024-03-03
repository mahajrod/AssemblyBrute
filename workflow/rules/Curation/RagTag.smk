

rule ragtag: #
    input:
        fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/scaffolds/{genome_prefix}.input.{haplotype}.fasta",
        reference_fasta=out_dir_path / "data/reference/{reference}/{reference}.softmasked.fasta"
    output:
        ragtag_fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/scaffolds/ragtag/{reference}/{genome_prefix}.{haplotype}.to.{reference}.fasta",
        ragtag_agp=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/scaffolds/ragtag/{reference}/{genome_prefix}.{haplotype}.to.{reference}.agp",
        ragtag_stats=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/scaffolds/ragtag/{reference}/{genome_prefix}.{haplotype}.to.{reference}.stats",
    log:
        ragtag=output_dict["log"]  / "create_gap_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{reference}.ragtag.log",
        ln=output_dict["log"]  / "create_gap_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{reference}.ln.log",
        cluster_log=output_dict["cluster_log"] / "create_gap_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{reference}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_gap_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{reference}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_gap_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{reference}.benchmark.txt"
    conda:
        config["conda"]["ragtag"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["ragtag"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("ragtag"),
        cpus=parameters["threads"]["ragtag"],
        time=parameters["time"]["ragtag"],
        mem=parameters["memory_mb"]["ragtag"]
    threads: parameters["threads"]["ragtag"]

    shell:
        " RAGTAG_DIR=`dirname {output.ragtag_fasta}`; "
        " ragtag.py scaffold -t {threads} -o ${{RAGTAG_DIR}}  -w {input.reference_fasta} {input.fasta} > {log.ragtag} 2>&1; "
        " ln -sf ragtag.scaffold.fasta {output.ragtag_fasta} > {log.ln} 2>&1; "
        " ln -sf ragtag.scaffold.agp {output.ragtag_agp} > {log.ln} 2>&1; "
        " ln -sf ragtag.scaffold.stats {output.ragtag_stats} > {log.ln} 2>&1; "