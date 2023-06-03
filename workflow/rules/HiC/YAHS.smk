
rule yahs: #
    input:
        bam=lambda wildcards: out_dir_path / "{0}/{1}/{2}/alignment/{3}/{4}.{5}.{3}.{2}.bwa.filtered.rmdup.bam".format(stage_dict["hic_scaffolding"]["prev_stage"],
                                                                                    wildcards.prev_stage_parameters, wildcards.haplotype,
                                                                                    stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..yahs_" + wildcards.hic_scaffolding_parameters]["option_set"]["phasing_kmer_length"],
                                                                                    wildcards.genome_prefix,
                                                                                    stage_dict["hic_scaffolding"]["prev_stage"]),
        reference=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta" % (stage_dict["hic_scaffolding"]["prev_stage"],
                                                                                                       stage_dict["hic_scaffolding"]["prev_stage"])) ,
        reference_fai=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta.fai" % (stage_dict["hic_scaffolding"]["prev_stage"],
                                                                                                         stage_dict["hic_scaffolding"]["prev_stage"]))
    output:
        fasta=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}_scaffolds_final.fa",
        alias=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{genome_prefix}.hic_scaffolding.{haplotype, [^.]+}.fasta",
    params:
        output_prefix=lambda wildcards: out_dir_path / "hic_scaffolding/{0}..yahs_{1}/{2}/scaffolding/{3}".format(wildcards.prev_stage_parameters,
                                                                                                             wildcards.hic_scaffolding_parameters,
                                                                                                             wildcards.haplotype,
                                                                                                             wildcards.genome_prefix),
        min_contig_len=lambda wildcards: stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..yahs_" + wildcards.hic_scaffolding_parameters]["option_set"]["min_contig_len"],
        min_mapping_quality=lambda wildcards: stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..yahs_" + wildcards.hic_scaffolding_parameters]["option_set"]["min_mapping_quality"],
        restriction_seq=parameters["tool_options"]["yahs"]["restriction_seq"][config["hic_enzyme_set"]],
    log:
        yahs=output_dict["log"]  / "yahs.hic_scaffolding.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.yahs.log",
        ln=output_dict["log"]  / "yahs.hic_scaffolding.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "yahs.hic_scaffolding.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "yahs.hic_scaffolding.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "yahs.hic_scaffolding.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["yahs"] ,
        time=parameters["time"]["yahs"],
        mem=parameters["memory_mb"]["yahs"]
    threads: parameters["threads"]["yahs"]

    shell:
        " yahs -q {params.min_mapping_quality} -l {params.min_contig_len} -o {params.output_prefix} "
        " {input.reference} {input.bam} > {log.yahs} 2>&1;"
        " ln {output.fasta} {output.alias} > {log.ln} 2>&1"
