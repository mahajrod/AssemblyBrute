
rule ragtag: #
    input:
        fasta=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta" % (stage_dict["ref_scaffolding"]["prev_stage"],
                                                                                                   stage_dict["ref_scaffolding"]["prev_stage"])),
        reference_fasta=out_dir_path / "data/reference/{reference}/{reference}.fasta",
        reference_syn=out_dir_path / "data/reference/{reference}/{reference}.syn",
    output:
        ragtag_fasta=out_dir_path / "ref_scaffolding/{prev_stage_parameters, [^/]+}..ragtag_{ref_scaf_parameters, [^/]+}@{reference, [^/]+}/{haplotype, [^/]+}/{genome_prefix, [^/]+}.ragtag.{haplotype}.fasta",
        ragtag_agp=out_dir_path / "ref_scaffolding/{prev_stage_parameters, [^/]+}..ragtag_{ref_scaf_parameters, [^/]+}@{reference, [^/]+}/{haplotype, [^/]+}/{genome_prefix, [^/]+}.ragtag.{haplotype}.agp",
        ragtag_stats=out_dir_path / "ref_scaffolding/{prev_stage_parameters, [^/]+}..ragtag_{ref_scaf_parameters, [^/]+}@{reference, [^/]+}/{haplotype, [^/]+}/{genome_prefix, [^/]+}.ragtag.{haplotype}.stats",
    params:
        min_aln_len=lambda wildcards: parse_option("min_aln_len", parameters["tool_options"]["ragtag"][wildcards.ref_scaf_parameters], " -f ")

    log:
        ragtag=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag.{genome_prefix}.ragtag.{haplotype}.log",
        ln=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag.{genome_prefix}.ragtag.{haplotype}.ln.log",
        awk=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/ragtag.{genome_prefix}.ragtag.{haplotype}.awk.log",
        cluster_log=output_dict["cluster_log"] / "ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.benchmark.txt"
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
        " ragtag.py scaffold -t {threads} -o ${{RAGTAG_DIR}} {params.min_aln_len} -w {input.reference_fasta} {input.fasta} > {log.ragtag} 2>&1; "
        " ln -sf ragtag.scaffold.fasta {output.ragtag_fasta} > {log.ln} 2>&1; "
        " ln -sf ragtag.scaffold.agp {output.ragtag_agp} >> {log.ln} 2>&1; "
        " ln -sf ragtag.scaffold.stats {output.ragtag_stats} >> {log.ln} 2>&1; "

rule rename_ragtag_scaffolds:
    input:
        fasta=rules.ragtag.output.ragtag_fasta,
        agp=rules.ragtag.output.ragtag_agp,
        stats=rules.ragtag.output.ragtag_stats,
        reference_syn=out_dir_path / "data/reference/{reference}/{reference}.syn",
    output:
        final_fasta=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference, [^/]+}/{genome_prefix, [^/]+}.ref_scaffolding.{haplotype, [^/]+}.fasta",
        final_agp=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference, [^/]+}/{genome_prefix, [^/]+}.ref_scaffolding.{haplotype, [^/]+}.agp",
        final_stats=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference, [^/]+}/{genome_prefix, [^/]+}.ref_scaffolding.{haplotype, [^/]+}.stats",
        ragtag_syn=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference, [^/]+}/{haplotype, [^/]+}/{genome_prefix}.ref_scaffolding.{haplotype}.syn",
    log:
        ln=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/rename_ragtag_scaffolds.{genome_prefix}.ref_scaffolding.{haplotype}.ln.log",
        awk=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/rename_ragtag_scaffolds.{genome_prefix}.ref_scaffolding.{haplotype}.awk.log",
        rename_fasta=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/rename_ragtag_scaffolds.{genome_prefix}.ref_scaffolding.{haplotype}.rename_fasta.log",
        rename_agp=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}/{haplotype}/rename_ragtag_scaffolds.{genome_prefix}.ref_scaffolding.{haplotype}.rename_agp.log",
        cluster_log=output_dict["cluster_log"] / "rename_ragtag_scaffolds.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "rename_ragtag_scaffolds{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "rename_ragtag_scaffolds.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}@{reference}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("rename_ragtag_scaffolds"),
        cpus=parameters["threads"]["rename_ragtag_scaffolds"],
        time=parameters["time"]["rename_ragtag_scaffolds"],
        mem=parameters["memory_mb"]["rename_ragtag_scaffolds"]
    threads: parameters["threads"]["rename_ragtag_scaffolds"]

    shell:
        " RAGTAG_DIR=`dirname {input.fasta}`; "
        " awk '{{print $2\"_RagTag\tps\"$1 }}' {input.reference_syn} > {output.ragtag_syn} 2>{log.awk}; "
        " rename_sequence_ids.py -i {input.fasta} -o {output.final_fasta} -l -s {output.ragtag_syn} -k 0 -c 1 > {log.rename_fasta} 2>&1; "
        " ln -sf {wildcards.haplotype}/{wildcards.genome_prefix}.ref_scaffolding.{wildcards.haplotype}.stats {output.final_stats} > {log.ln} 2>&1; "
        " replace_column_value_by_syn.py -i {input.agp} -o {output.final_agp} -s {output.ragtag_syn} -c 0 > {log.rename_agp} 2>&1; " # -k 0 -c 1