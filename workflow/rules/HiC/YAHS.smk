
rule yahs: #
    input:
        hic_qc=lambda wildcards: expand(out_dir_path / "{0}/{1}/{2}/alignment/{3}/{4}.{0}.{3}.{2}.rmdup.inter_{5}.txt".format(stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..yahs_" + wildcards.hic_scaffolding_parameters]["prev_stage"],
                                                                                                                              wildcards.prev_stage_parameters,
                                                                                                                              wildcards.haplotype,
                                                                                                                              stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..yahs_" + wildcards.hic_scaffolding_parameters]["option_set"]["phasing_kmer_length"],
                                                                                                                              wildcards.genome_prefix, "{mapq}"),
                                        mapq=parameters["tool_options"]["juicer_tools_qc"]["mapping_quality"],
                                        allow_missing=True),
        bam=lambda wildcards: out_dir_path / "{0}/{1}/{2}/alignment/{3}/{4}.{0}.{3}.{2}.rmdup.bam".format(stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..yahs_" + wildcards.hic_scaffolding_parameters]["prev_stage"],#stage_dict["hic_scaffolding"]["prev_stage"],
                                                                                                          wildcards.prev_stage_parameters,
                                                                                                          wildcards.haplotype,
                                                                                                          stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..yahs_" + wildcards.hic_scaffolding_parameters]["option_set"]["phasing_kmer_length"],
                                                                                                          wildcards.genome_prefix),
        reference=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta" % (stage_dict["hic_scaffolding"]["prev_stage"],
                                                                                                       stage_dict["hic_scaffolding"]["prev_stage"])) ,
        reference_fai=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta.fai" % (stage_dict["hic_scaffolding"]["prev_stage"],
                                                                                                               stage_dict["hic_scaffolding"]["prev_stage"]))
    output:
        fasta=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}_scaffolds_final.fa",
        bin=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.bin",
        agp=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}_scaffolds_final.agp",
        alias=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{genome_prefix}.hic_scaffolding.{haplotype, [^.]+}.fasta",
    params:
        #output_prefix=lambda wildcards: out_dir_path / "hic_scaffolding/{0}..yahs_{1}/{2}/scaffolding/{3}".format(wildcards.prev_stage_parameters,
        #                                                                                                     wildcards.hic_scaffolding_parameters,
        #                                                                                                     wildcards.haplotype,
        #                                                                                                     wildcards.genome_prefix),
        min_contig_len=lambda wildcards: stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..yahs_" + wildcards.hic_scaffolding_parameters]["option_set"]["min_contig_len"],
        min_mapping_quality=lambda wildcards: stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..yahs_" + wildcards.hic_scaffolding_parameters]["option_set"]["min_mapping_quality"],
        #restriction_seq=parameters["tool_options"]["yahs"]["restriction_seq"][config["hic_enzyme_set"]],
        restriction_seq=parse_option(config["hic_enzyme_set"], parameters["tool_options"]["yahs"]["restriction_seq"], "-e",)
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
        " OUTPUT_PREFIX={output.fasta}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%_scaffolds_final.fa}}; "
        " yahs -q {params.min_mapping_quality} -l {params.min_contig_len} -o ${{OUTPUT_PREFIX}} "
        " {input.reference} {input.bam} > {log.yahs} 2>&1;"
        " ln -sf {wildcards.haplotype}/scaffolding/`basename {output.fasta}` {output.alias} > {log.ln} 2>&1"

rule yahs_juicer_pre: #
    input:
        reference_fai=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta.fai" % (stage_dict["hic_scaffolding"]["prev_stage"],
                                                                                                               stage_dict["hic_scaffolding"]["prev_stage"])),
        bin=rules.yahs.output.bin,
        agp=rules.yahs.output.agp
    output:
        links_bed=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.bed",
        liftover_agp=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.liftover.agp",
        assembly=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.assembly",
        assembly_agp=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.assembly.agp",
        log=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.log",
    log:
        #std=output_dict["log"]  / "yahs_juicer_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.log",
        mv=output_dict["log"]  / "yahs_juicer_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.mv.log",
        cluster_log=output_dict["cluster_log"] / "yahs_juicer_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "yahs_juicer_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "yahs_juicer_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["yahs_juicer_pre"] ,
        time=parameters["time"]["yahs_juicer_pre"],
        mem=parameters["memory_mb"]["yahs_juicer_pre"]
    threads: parameters["threads"]["yahs_juicer_pre"]

    shell:
        " OUTPUT_PREFIX={output.links_bed}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.bed}}; "
        " juicer pre -a -o ${{OUTPUT_PREFIX}} {input.bin} {input.agp} {input.reference_fai} > {output.log} 2>&1;"
        " mv ${{OUTPUT_PREFIX}}.txt {output.links_bed} > {log.mv} 2>&1; "

rule juicer_tools_pre: #
    input:
        yahs_juicer_pre_log=rules.yahs_juicer_pre.output.log,
        yahs_juicer_pre_bed=rules.yahs_juicer_pre.output.links_bed
    output:
        hic=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.hic",
    log:
        juicer=output_dict["log"]  / "juicer_tools_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.juicer.log",
        cat=output_dict["log"]  / "juicer_tools_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cat.log",
        grep=output_dict["log"]  / "juicer_tools_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.grep.log",
        awk=output_dict["log"]  / "juicer_tools_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.awk.log",
        cluster_log=output_dict["cluster_log"] / "juicer_tools_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "juicer_tools_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "juicer_tools_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["juicer_tools_pre"] ,
        time=parameters["time"]["juicer_tools_pre"],
        mem=parameters["memory_mb"]["juicer_tools_pre"]
    threads: parameters["threads"]["juicer_tools_pre"]

    shell: # juicer_tools elder than 1.9.9 seems to be incompartible with yahs
        " java -jar -Xmx{resources.mem}m workflow/external_tools/juicer/juicer_tools.1.9.9_jcuda.0.8.jar pre " #--threads {threads}  
        " {input.yahs_juicer_pre_bed} {output.hic} <(cat {input.yahs_juicer_pre_log} 2>{log.cat} | "
        " grep PRE_C_SIZE 2>{log.grep} | awk '{{print $2\" \"$3}}' 2>{log.awk}) > {log.juicer} 2>&1; "
