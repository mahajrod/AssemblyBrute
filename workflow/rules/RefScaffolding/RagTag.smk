

rule ragtag: #
    input:
        fasta=out_dir_path / ("{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta" % stage_dict["ref_scaffolding"]["prev_stage"]),
        reference_fasta=out_dir_path / "data/reference/{reference}/{reference}.softmasked.fasta"
    output:
        alias_ragtag_fasta=out_dir_path / "ref_scaffolding/{prev_stage_parameters, [^/]}..ragtag_{ref_scaf_parameters, [^/]}_{reference, [^/]}/{genome_prefix, [^/]}.ref_scaffolding.{haplotype, [^/]}.fasta",
        alias_ragtag_agp=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}_{reference}/{genome_prefix}.ref_scaffolding.{haplotype}.agp",
        alias_ragtag_stats=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}_{reference}/{genome_prefix}.ref_scaffolding.{haplotype}.stats",
        ragtag_fasta=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}_{reference}/{haplotype, [^.]+}/{genome_prefix, [^/]}.ref_scaffolding.{haplotype, [^/]}.fasta",
        ragtag_agp=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}_{reference}/{haplotype, [^.]+}/{genome_prefix, [^/]}.ref_scaffolding.{haplotype, [^/]}.agp",
        ragtag_stats=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}_{reference}/{haplotype, [^.]+}/{genome_prefix, [^/]}.ref_scaffolding.{haplotype, [^/]}.stats",
    log:
        ragtag=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}_{reference}/{haplotype}/ragtag.{genome_prefix}.ref_scaffolding.{haplotype}.ragtag.log",
        ln=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}_{reference}/{haplotype}/ragtag.{genome_prefix}.ref_scaffolding.{haplotype}.ragtag.ln.log",
        cluster_log=output_dict["cluster_log"] / "ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}_{reference}.{genome_prefix}.{haplotype}.{reference}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}_{reference}.{genome_prefix}.{haplotype}.{reference}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}_{reference}.{genome_prefix}.{haplotype}.benchmark.txt"
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
        " ln -sf ragtag.scaffold.agp {output.ragtag_agp} >> {log.ln} 2>&1; "
        " ln -sf ragtag.scaffold.stats {output.ragtag_stats} >> {log.ln} 2>&1; "
        " ln -sf {wildcards.haplotype}/{wildcards.genome_prefix}.ref_scaffolding.{wildcards.haplotype}.fasta {output.alias_ragtag_fasta} >> {log.ln} 2>&1; "
        " ln -sf {wildcards.haplotype}/{wildcards.genome_prefix}.ref_scaffolding.{wildcards.haplotype}.agp {output.alias_ragtag_agp} >> {log.ln} 2>&1; "
        " ln -sf {wildcards.haplotype}/{wildcards.genome_prefix}.ref_scaffolding.{wildcards.haplotype}.stats {output.alias_ragtag_stats} >> {log.ln} 2>&1; "
