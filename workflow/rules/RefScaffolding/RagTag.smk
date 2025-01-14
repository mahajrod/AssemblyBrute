

def get_reference_id_from_ref_scaffolding_parameters(parameters):
    return parameters.split("!")[-1]

def get_reference_fasta_path_from_ref_scaffolding_parameters(parameters):
    ref_id = get_reference_id_from_ref_scaffolding_parameters(parameters)
    return out_dir_path / ("data/reference/{0}/{0}.softmasked.fasta".format(ref_id))

rule ragtag: #
    input:
        fasta=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta" %(stage_dict["ref_scaffolding"]["prev_stage"],
                                                                                                  stage_dict["ref_scaffolding"]["prev_stage"])),
        reference_fasta=lambda wildcards: get_reference_fasta_path_from_ref_scaffolding_parameters(wildcards.ref_scaf_parameters) #out_dir_path / "data/reference/{reference}/{reference}.softmasked.fasta"
    output:
        #fasta_prefix = results / ref_scaffolding / hifiasm_l3primary_no_hic..purge_dups_keep_hicov..ragtag_default!bCatUst1.pri.v2 / TurFus1.ref_scaffolding.hap0
        alias_ragtag_fasta=out_dir_path / "ref_scaffolding/{prev_stage_parameters, [^/]}..ragtag_{ref_scaf_parameters, [^/]}/{genome_prefix, [^/]}.ref_scaffolding.{haplotype, [^/]}.fasta",
        alias_ragtag_agp=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}/{genome_prefix}.ref_scaffolding.{haplotype}.agp",
        alias_ragtag_stats=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}/{genome_prefix}.ref_scaffolding.{haplotype}.stats",
        ragtag_fasta=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}/{haplotype, [^.]+}/{genome_prefix, [^/]}.ref_scaffolding.{haplotype}.fasta",
        ragtag_agp=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}/{haplotype, [^.]+}/{genome_prefix, [^/]}.ref_scaffolding.{haplotype}.agp",
        ragtag_stats=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}/{haplotype, [^.]+}/{genome_prefix, [^/]}.ref_scaffolding.{haplotype}.stats",
    log:
        ragtag=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}/{haplotype}/ragtag.{genome_prefix}.ref_scaffolding.{haplotype}.ragtag.log",
        ln=out_dir_path / "ref_scaffolding/{prev_stage_parameters}..ragtag_{ref_scaf_parameters}/{haplotype}/ragtag.{genome_prefix}.ref_scaffolding.{haplotype}.ragtag.ln.log",
        cluster_log=output_dict["cluster_log"] / "ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "ragtag.{prev_stage_parameters}..ragtag_{ref_scaf_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
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
