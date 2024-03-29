#OBSOLETE CODE

rule salsa2: #
    input:
        bed=lambda wildcards: out_dir_path / "{0}/{1}/{2}/alignment/{3}/{4}.{5}.{3}.{2}.rmdup.bed".format(stage_dict["hic_scaffolding"]["prev_stage"],
                                                                                    wildcards.prev_stage_parameters, wildcards.haplotype,
                                                                                    stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..salsa2_" + wildcards.hic_scaffolding_parameters]["option_set"]["phasing_kmer_length"],
                                                                                    wildcards.genome_prefix,
                                                                                    stage_dict["hic_scaffolding"]["prev_stage"]),
        reference=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta" % (stage_dict["hic_scaffolding"]["prev_stage"],
                                                                                                       stage_dict["hic_scaffolding"]["prev_stage"])) ,
        reference_fai=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta.fai" % (stage_dict["hic_scaffolding"]["prev_stage"],
                                                                                                         stage_dict["hic_scaffolding"]["prev_stage"]))
        #reference=out_dir_path  / ("purge_dups/{assembler}/%s.purge_dups.{assembler}.{haplotype}.fasta" % config["genome_name"]),
        #reference_fai=out_dir_path  / ("purge_dups/{assembler}/%s.purge_dups.{assembler}.{haplotype}.fai" % config["genome_name"])
    output:
        dir=directory(out_dir_path / "hic_scaffolding/{prev_stage_parameters}..salsa2_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}"),
        #dir=directory(out_dir_path  / ("hic_scaffolding/{assembler}/{haplotype}/scaffolding/%s/"  % config["genome_name"])),
        fasta=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..salsa2_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}/scaffolds_FINAL.fasta",
        #fasta=out_dir_path  / ("hic_scaffolding/{assembler}/{haplotype}/scaffolding/%s/scaffolds_FINAL.fasta"  % config["genome_name"]),
        alias=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..salsa2_{hic_scaffolding_parameters}/{genome_prefix}.hic_scaffolding.{haplotype, [^.]+}.fasta",
        #alias=out_dir_path  / ("hic_scaffolding/{assembler}/%s.hic_scaffolding.{assembler}.{haplotype}.fasta" % config["genome_name"])
    params:
        min_contig_len=lambda wildcards: stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..salsa2_" + wildcards.hic_scaffolding_parameters]["option_set"]["min_contig_len"],
        restriction_seq=parameters["tool_options"]["salsa2"]["restriction_seq"][config["hic_enzyme_set"]],
    log:
        salsa=output_dict["log"]  / "salsa2.hic_scaffolding.{prev_stage_parameters}..salsa2_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.salsa.log",
        ln=output_dict["log"]  / "salsa2.hic_scaffolding.{prev_stage_parameters}..salsa2_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "salsa2.hic_scaffolding.{prev_stage_parameters}..salsa2_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "salsa2.hic_scaffolding.{prev_stage_parameters}..salsa2_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "salsa2.hic_scaffolding.{prev_stage_parameters}..salsa2_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["salsa2"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["salsa2"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("salsa2"),
        cpus=parameters["threads"]["salsa2"] ,
        time=parameters["time"]["salsa2"],
        mem=parameters["memory_mb"]["salsa2"]
    threads: parameters["threads"]["salsa2"]

    shell:
        " run_pipeline.py -a {input.reference} -l {input.reference_fai} -b {input.bed} -e {params.restriction_seq} "
        " -m yes -p yes "
        " -o {output.dir} > {log.salsa} 2>&1;"
        " ln -sf {wildcards.haplotype}/scaffolding/{wildcards.genome_prefix}/`basename {output.fasta}` {output.alias} > {log.ln} 2>&1"
