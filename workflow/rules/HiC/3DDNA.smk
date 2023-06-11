
rule generate_site_positions: #
    input:
        fasta=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..3ddna_" + wildcards.hic_scaffolding_parameters]["prev_stage"],
                                                                                      wildcards.prev_stage_parameters,
                                                                                      wildcards.genome_prefix,
                                                                                      wildcards.haplotype)
    output:
        #alias_fasta=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.input.{haplotype}.fasta",
        restriction_site_file=out_dir_path / ("hic_scaffolding/{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.input.{haplotype}_%s.txt" % ["hic_enzyme_set"]),
    params:
        restriction_seq=config["hic_enzyme_set"]
    log:
        #ln=output_dict["log"]  / "generate_site_positions.hic_scaffolding.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.ln.log",
        #ln=output_dict["log"]  / "generate_site_positions.hic_scaffolding.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.ln.log",
        sites=output_dict["log"]  / "generate_site_positions.hic_scaffolding.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.sites.log",
        cluster_log=output_dict["cluster_log"] / "generate_site_positions.hic_scaffolding.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "generate_site_positions.hic_scaffolding.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "generate_site_positions.hic_scaffolding.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["generate_site_positions"] ,
        time=parameters["time"]["generate_site_positions"],
        mem=parameters["memory_mb"]["generate_site_positions"]
    threads: parameters["threads"]["generate_site_positions"]

    shell:
        #" ln -sf ../../../../../../{input.fasta} {output.alias_fasta} > {log.ln} 2>&1; "
        " OUTPUT_PREFIX={output.restriction_site_file}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%_{params.restriction_seq}.txt}; "
        " ./workflow/external_tools/juicer/misc/generate_site_positions.py {params.restriction_seq} ${{OUTPUT_PREFIX}} {input.fasta} > {log.sites} 2>&1; "

rule juicer: #
    input:
        fasta=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..3ddna_" + wildcards.hic_scaffolding_parameters]["prev_stage"],
                                                                                  wildcards.prev_stage_parameters,
                                                                                  wildcards.genome_prefix,
                                                                                  wildcards.haplotype),
        index=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta.bwt".format(stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..3ddna_" + wildcards.hic_scaffolding_parameters]["prev_stage"],
                                                                                      wildcards.prev_stage_parameters,
                                                                                      wildcards.genome_prefix,
                                                                                      wildcards.haplotype),
        restriction_site_file=rules.generate_site_positions.output.restriction_site_file if config["hic_enzyme_set"] != "OmniC" else [],
        fastqs=lambda wildcards: expand(output_dict["data"] / "fastq/hic/raw/{0}{1}".format("{fileprefix}", config["fastq_extension"]) if parameters["tool_options"]["3ddna"][wildcards.prev_stage_parameters + "..3ddna_" + wildcards.hic_scaffolding_parameters]["phasing_kmer_length"] == "NA" else \
                                        out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}".format(config["phasing_stage"], #wildcards.assembly_stage,
                                                                                                 detect_phasing_parameters(wildcards.prev_stage_parameters + "..3ddna_" + wildcards.hic_scaffolding_parameters,
                                                                                                                           config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
                                                                                                 wildcards.haplotype,
                                                                                                 parameters["tool_options"]["3ddna"][wildcards.prev_stage_parameters + "..3ddna_" + wildcards.hic_scaffolding_parameters]["phasing_kmer_length"],
                                                                                                 "{fileprefix}",
                                                                                                 config["fastq_extension"]
                                                                                                 ),
                                        fileprefix=input_file_prefix_dict["hic"])
    params:
        restriction_seq=config["hic_enzyme_set"] if config["hic_enzyme_set"] != "OmniC" else "none"
    output:
        merged_no_dups=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}_scaffolds_final.fa",
    log:
        juicer=output_dict["log"]  / "juicer.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.juicer.log",
        mkdir=output_dict["log"]  / "juicer.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.mkdir.log",
        ln=output_dict["log"]  / "juicer.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "juicer.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "juicer.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "juicer.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["juicer"] ,
        time=parameters["time"]["juicer"],
        mem=parameters["memory_mb"]["juicer"]
    threads: parameters["threads"]["juicer"]

    shell:
        " OUTPUT_DIR"
        " mkdir -p ${{OUTPUT_DIR}}/fastq > {log.mkdir} 2>&1 ;"
        " > {log.ln}; "
        " for FILE in {input.fastqs};"
        " do"
        "   ln -sf ../../../../../../${{FILE}} ${{OUTPUT_DIR}}/fastq >> {log.ln} 2>&1; "
        " done; "
        " cd ${{OUTPU_DIR}}; "
        " ../../../../../workflow/external_tools/juicer/scripts/juicer.sh -t {threads} -D `realpath ../../../../../workflow/external_tools/juicer/` "
        " -g {wildcards.genome_prefix} -s {params.restriction_seq} –z `realpath ../../../../../{input.fasta}` –y `basename {input.restriction_site_file}` "
        " –p assembly -d `realpath ./` > ../../../../../{log.juicer} 2>&1; "

"""
rule juicer_tools_pre: #
    input:
        yahs_juicer_pre_log=rules.yahs_juicer_pre.output.log,
        yahs_juicer_pre_bed=rules.yahs_juicer_pre.output.links_bed
    output:
        hic=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.hic",
    log:
        juicer=output_dict["log"]  / "juicer_tools_pre.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.juicer.log",
        cat=output_dict["log"]  / "juicer_tools_pre.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cat.log",
        grep=output_dict["log"]  / "juicer_tools_pre.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.grep.log",
        awk=output_dict["log"]  / "juicer_tools_pre.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.awk.log",
        cluster_log=output_dict["cluster_log"] / "juicer_tools_pre.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "juicer_tools_pre.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "juicer_tools_pre.{prev_stage_parameters}..3ddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["juicer_tools_pre"] ,
        time=parameters["time"]["juicer_tools_pre"],
        mem=parameters["memory_mb"]["juicer_tools_pre"]
    threads: parameters["threads"]["juicer_tools_pre"]

    shell: # juicer_tools elder than 1.9.9 seems to be incompartible with yahs
        " java -jar -Xmx{resources.mem}m workflow/external_tools/juicer/juicer_tools.jar pre " #--threads {threads}  
        " {input.yahs_juicer_pre_bed} {output.hic} <(cat {input.yahs_juicer_pre_log} 2>{log.cat} | "
        " grep PRE_C_SIZE 2>{log.grep} | awk '{{print $2\" \"$3}}' 2>{log.awk}) > {log.juicer} 2>&1; "
"""