"""
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
"""

rule juicer: #
    input:
        fasta=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["prev_stage"],
                                                                                  wildcards.prev_stage_parameters,
                                                                                  wildcards.genome_prefix,
                                                                                  wildcards.haplotype),
        index=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta.bwt".format(stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["prev_stage"],
                                                                                      wildcards.prev_stage_parameters,
                                                                                      wildcards.genome_prefix,
                                                                                      wildcards.haplotype),
        restriction_site_file=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}_{4}.txt".format(stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["prev_stage"],
                                                                                                    wildcards.prev_stage_parameters,
                                                                                                    wildcards.genome_prefix,
                                                                                                    wildcards.haplotype,
                                                                                                    config[["hic_enzyme_set"]]) if config["hic_enzyme_set"] not in config["no_motif_enzyme_sets"] else [],
        forward_fastqs=lambda wildcards: expand(output_dict["data"] / "fastq/hic/raw/{0}{1}".format("{fileprefix}", config["fastq_extension"]) if parameters["tool_options"]["threeddna"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["phasing_kmer_length"] == "NA" else \
                                        out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}".format(config["phasing_stage"], #wildcards.assembly_stage,
                                                                                                 detect_phasing_parameters(wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters,
                                                                                                                           config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
                                                                                                 wildcards.haplotype,
                                                                                                 parameters["tool_options"]["threeddna"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["phasing_kmer_length"],
                                                                                                 "{fileprefix}",
                                                                                                 config["fastq_extension"]
                                                                                                 ),
                                        fileprefix=input_file_prefix_dict["hic"][::2]),
        reverse_fastqs=lambda wildcards: expand(output_dict["data"] / "fastq/hic/raw/{0}{1}".format("{fileprefix}", config["fastq_extension"]) if parameters["tool_options"]["threeddna"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["phasing_kmer_length"] == "NA" else \
                                        out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}".format(config["phasing_stage"], #wildcards.assembly_stage,
                                                                                                 detect_phasing_parameters(wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters,
                                                                                                                           config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
                                                                                                 wildcards.haplotype,
                                                                                                 parameters["tool_options"]["threeddna"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["phasing_kmer_length"],
                                                                                                 "{fileprefix}",
                                                                                                 config["fastq_extension"]
                                                                                                 ),
                                        fileprefix=input_file_prefix_dict["hic"][1::2])
    params:
        restriction_seq=config["hic_enzyme_set"]  if config["hic_enzyme_set"] not in config["no_motif_enzyme_sets"] else "none",
        fastq_extensions=config["fastq_extension"]
    output:
        merged_no_dups=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.merged_nodups.txt",
        merged_dedup_bam=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.merged_dedup.bam",
        merged_inter_30=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.inter_30.txt",
        merged_inter=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.inter.txt",
    log:
        juicer=output_dict["log"]  / "juicer.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.juicer.log",
        mkdir=output_dict["log"]  / "juicer.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.mkdir.log",
        ln=output_dict["log"]  / "juicer.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "juicer.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "juicer.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "juicer.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["juicer"] ,
        time=parameters["time"]["juicer"],
        mem=parameters["memory_mb"]["juicer"]
    threads: parameters["threads"]["juicer"]

    shell:
        " OUTPUT_DIR=`dirname {output.merged_no_dups}; `"
        " mkdir -p ${{OUTPUT_DIR}}/fastq > {log.mkdir} 2>&1 ; "
        " > {log.ln}; "
        " for FILE in {input.forward_fastqs}; "
        " do"
        "       BASE_FILENAME=`basename ${{FILE}}`; "    
        "       ln ${{FILE}} ${{OUTPUT_DIR}}/fastq/${{BASE_FILENAME%{params.fastq_extensions}}}_R1_001{params.fastq_extensions} >> {log.ln} 2 >&1; "
        " done;"
        " for FILE in {input.reverse_fastqs}; "
        " do"
        "       BASE_FILENAME=`basename ${{FILE}}`; "    
        "       ln ${{FILE}} ${{OUTPUT_DIR}}/fastq/${{BASE_FILENAME%{params.fastq_extensions}}}_R2_001{params.fastq_extensions} >> {log.ln} 2 >&1; "
        " done;"
        " SCRIPT=`realpath ./workflow/external_tools/juicer/scripts/juicer.sh`;"
        " JUICER_DIR=`realpath ./workflow/external_tools/juicer/`; "
        " FASTA=`realpath {input.fasta}`; "
        " RESTRICTION_SITE_FILE=`realpath {input.restriction_site_file}`; "
        " JUICER_LOG=`realpath {log.juicer}`; "
        " ${{SCRIPT}} -t {threads} -D ${{JUICER_DIR}} -g {wildcards.genome_prefix} -s {params.restriction_seq} "
        " –z ${{FASTA}} –y ${{RESTRICTION_SITE_FILE}} –p assembly -d ${{OUTPUT_DIR}} > ${{JUICER_LOG}} 2>&1; "
        " mv ${{OUTPUT_DIR}}/fastq/merged_nodups.txt {output.merged_no_dups}; "
        " mv ${{OUTPUT_DIR}}/fastq/merged_dedup.bam {output.merged_dedup_bam}; " 
        " mv ${{OUTPUT_DIR}}/fastq/inter_30.txt {output.merged_inter_30}; "
        " mv ${{OUTPUT_DIR}}/fastq/inter.txt {output.merged_inter}; " 
        " rm -r ${{OUTPUT_DIR}}/splits ${{OUTPUT_DIR}}/aligned; "

rule threeddna: #
    input:
        fasta=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["prev_stage"],
                                                                                  wildcards.prev_stage_parameters,
                                                                                  wildcards.genome_prefix,
                                                                                  wildcards.haplotype),
        merged_nodups=rules.juicer.output.merged_no_dups
    params:
        restriction_seq=config["hic_enzyme_set"]  if config["hic_enzyme_set"] not in config["no_motif_enzyme_sets"] else "none",
        fastq_extensions=config["fastq_extension"]
    output:
        draft_fasta=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.input.{haplotype}.fasta",
        rawchrom_hic=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.input.{haplotype}.rawchrom.hic",
        rawchrom_assembly=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.input.{haplotype}.rawchrom.assembly",
        hic_fasta=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.input.{haplotype}_HiC.fasta",
        alias_fasta=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}/{genome_prefix}.hic_scaffolding.{haplotype, [^.]+}.fasta",

    log:
        threeddna=output_dict["log"]  / "threeddna.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.threeddna.log",
        ln=output_dict["log"]  / "threeddna.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "threeddna.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "threeddna.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "threeddna.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["threeddna"] ,
        time=parameters["time"]["threeddna"],
        mem=parameters["memory_mb"]["threeddna"]
    threads: parameters["threads"]["threeddna"]

    shell:
        " OUTPUT_DIR=`dirname {output}; ` "
        " SCRIPT=`realpath workflow/external_tools/3d-dna/run-3ddna-pipeline.sh`; "
        " INPUT_FASTA=`realpath {input.fasta}`; "
        " THREEDDNA_LOG=`realpath {log.threeddna}`; "
        " LN_LOG=`realpath {log.ln}`; "
        " > ${{LN_LOG}}; "
        " cd ${{OUTPUT_DIR}}; " 
        " ln ${{INPUT_FASTA}} {output.draft_fasta} >> ${{LN_LOG}} 2>&1; "
        " ${{SCRIPT}} `basename {output.draft_fasta}` `basename {input.merged_nodups}` > ${{THREEDDNA_LOG}} 2>&1; "
        " ln {wildcards.haplotype}/scaffolding/{wildcards.genome_prefix}.input.{wildcards.haplotype}_HiC.fasta "
        " ../../`basename {output.alias_fasta}` >> ${{LN_LOG}} 2>&1"
