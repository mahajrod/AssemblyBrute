localrules: create_fastq_links_for_juicer
ruleorder: create_fastq_links_for_juicer > create_fastq_links
#stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["option_set"]

#input_pairprefix_dict["hic"]
#input_forward_suffix_dict["hic"]
#input_reverse_suffix_dict["hic"]
def get_hic_reads_for_juicer(wildcards):
    output_directory = out_dir_path / "hic_scaffolding/{0}..threeddna_{1}/{2}/scaffolding/fastq/".format(wildcards.prev_stage_parameters,
                                                                                                         wildcards.hic_scaffolding_parameters,
                                                                                                         wildcards.haplotype)
    output_forward_suffix = "_R1_001"
    output_reverse_suffix = "_R2_001"
    if stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["option_set"]["phasing_kmer_length"] == "NA":
        forward_suffix = input_forward_suffix_dict["hic"]
        reverse_suffix = input_reverse_suffix_dict["hic"]
        directory = output_dict["data"] / "fastq/hic/raw/"
    else:
        forward_suffix = "_1"
        reverse_suffix = "_2"
        directory = out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/".format(config["phasing_stage"],
                                                                      detect_phasing_parameters(wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters,
                                                                                                config["phasing_stage"], stage_separator=".."),
                                                                      wildcards.haplotype,
                                                                      stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["option_set"]["phasing_kmer_length"],
                                                                      )
    input_forward_filelist = []
    input_reverse_filelist = []
    output_forward_filelist = []
    output_reverse_filelist = []

    for pairprefix in input_pairprefix_dict["hic"]:
        input_forward_filelist.append("{0}/{1}{2}{3}".format(directory, pairprefix, forward_suffix,config["fastq_extension"]))
        input_reverse_filelist.append("{0}/{1}{2}{3}".format(directory, pairprefix, reverse_suffix,config["fastq_extension"]))
        output_forward_filelist.append("{0}/{1}{2}{3}".format(output_directory, pairprefix, output_forward_suffix,config["fastq_extension"]))
        output_reverse_filelist.append("{0}/{1}{2}{3}".format(output_directory, pairprefix, output_reverse_suffix,config["fastq_extension"]))

    return input_forward_filelist, input_reverse_filelist, output_forward_filelist, output_reverse_filelist

"""
def get_reverse_reads(wildcards):
    reverse_reads = expand(output_dict["data"] / "fastq/hic/raw/{0}{1}".format("{fileprefix}", config["fastq_extension"]) if stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["option_set"]["phasing_kmer_length"] == "NA" else \
                                        out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}".format(config["phasing_stage"], #wildcards.assembly_stage,
                                                                                                 detect_phasing_parameters(wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters,
                                                                                                                           config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
                                                                                                 wildcards.haplotype,
                                                                                                 stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["option_set"]["phasing_kmer_length"],
                                                                                                 "{fileprefix}",
                                                                                                 config["fastq_extension"]
                                                                                                 ),
                                        fileprefix=input_file_prefix_dict["hic"][1::2])

    return reverse_reads
"""
rule create_fastq_links_for_juicer:
    input:
        forward_fastqs=lambda wildcards: get_hic_reads_for_juicer(wildcards)[0],
        reverse_fastqs=lambda wildcards: get_hic_reads_for_juicer(wildcards)[1],

    output:
        fastq_dir=directory(out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..threeddna_{hic_scaffolding_parameters, [^/]+}/{haplotype, [^.]+}/scaffolding/fastq")
    log:
        ln=output_dict["log"]  / "create_fastq_links_for_juicer.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "create_fastq_links_for_juicer.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_fastq_links_for_juicer.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_fastq_links_for_juicer.{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}.{haplotype}.benchmark.txt"
    #conda:
    #    config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("create_fastq_links_for_juicer"),
        cpus=parameters["threads"]["create_fastq_links_for_juicer"] ,
        time=parameters["time"]["create_fastq_links_for_juicer"],
        mem=parameters["memory_mb"]["create_fastq_links_for_juicer"]
    threads: parameters["threads"]["create_fastq_links_for_juicer"]

    run:
        shell(" > {0}; ".format(log.ln))
        shell("mkdir -p {0} ; ".format(output.fastq_dir))
        for input_forward, input_reverse, output_forward, output_reverse in zip(*get_hic_reads_for_juicer(wildcards)):
            shell("ln -sf `realpath {0}` {1} >> {2} 2>&1; ".format(input_forward, output_forward, log.ln))
            shell("ln -sf `realpath {0}` {1} >> {2} 2>&1; ".format(input_reverse, output_reverse, log.ln))


rule juicer: #
    input:
        fasta=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["prev_stage"],
                                                                                  wildcards.prev_stage_parameters,
                                                                                  wildcards.genome_prefix,
                                                                                  wildcards.haplotype),
        index=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta.ann".format(stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["prev_stage"],
                                                                                      wildcards.prev_stage_parameters,
                                                                                      wildcards.genome_prefix,
                                                                                      wildcards.haplotype),
        restriction_site_file=lambda wildcards: (out_dir_path / "{0}/{1}/{2}.{0}.{3}_{4}.txt".format(stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["prev_stage"],
                                                                                                     wildcards.prev_stage_parameters,
                                                                                                     wildcards.genome_prefix,
                                                                                                     wildcards.haplotype,
                                                                                                     config["hic_enzyme_set"]) ) if config["hic_enzyme_set"] not in config["no_motif_enzyme_sets"] else [],
        fastq_dir=rules.create_fastq_links_for_juicer.output.fastq_dir
        #forward_fastqs=lambda wildcards: get_hic_reads(wildcards)[0],
        #reverse_fastqs=lambda wildcards: get_hic_reads(wildcards)[1],
    params:
        restriction_seq=config["hic_enzyme_set"]  if config["hic_enzyme_set"] not in config["no_motif_enzyme_sets"] else "none",
        fastq_extensions=config["fastq_extension"]
    output:
        merged_no_dups=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..threeddna_{hic_scaffolding_parameters, [^/]+}/{haplotype, [^.]+}/scaffolding/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype}.merged_nodups.txt",
        merged_dedup_bam=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..threeddna_{hic_scaffolding_parameters, [^/]+}/{haplotype, [^.]+}/scaffolding/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype}.merged_dedup.bam",
        merged_inter_30=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..threeddna_{hic_scaffolding_parameters, [^/]+}/{haplotype, [^.]+}/scaffolding/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype}.inter_30.txt",
        merged_inter=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..threeddna_{hic_scaffolding_parameters, [^/]+}/{haplotype, [^.]+}/scaffolding/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype}.inter.txt",
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
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("juicer"),
        cpus=parameters["threads"]["juicer"] ,
        time=parameters["time"]["juicer"],
        mem=parameters["memory_mb"]["juicer"]
    threads: parameters["threads"]["juicer"]

    shell:
        " OUTPUT_DIR=`dirname {output.merged_no_dups}`; "
        " OUTPUT_DIR=`realpath -s ${{OUTPUT_DIR}}`; "
        " > {log.ln}; "
        " SCRIPT=`realpath -s ./workflow/external_tools/juicer/scripts/juicer.sh`;"
        " JUICER_DIR=`realpath -s ./workflow/external_tools/juicer/`; "
        " FASTA=`realpath -s {input.fasta}`; "
        " RESTRICTION_SITE_FILE=`realpath -s {input.restriction_site_file}`; "
        " ${{SCRIPT}} -t {threads} -T {threads} -D ${{JUICER_DIR}} -g {wildcards.genome_prefix} -s {params.restriction_seq} "
        " -z ${{FASTA}} -y ${{RESTRICTION_SITE_FILE}} --assembly -d ${{OUTPUT_DIR}} > {log.juicer} 2>&1; "
        " mv ${{OUTPUT_DIR}}/aligned/merged_nodups.txt {output.merged_no_dups}; "
        " mv ${{OUTPUT_DIR}}/aligned/merged_dedup.bam {output.merged_dedup_bam}; " 
        " mv ${{OUTPUT_DIR}}/aligned/inter_30.txt {output.merged_inter_30}; "
        " mv ${{OUTPUT_DIR}}/aligned/inter.txt {output.merged_inter}; " 
        " rm -r ${{OUTPUT_DIR}}/splits ; " # ${{OUTPUT_DIR}}/aligned
        #
        #" mkdir -p ${{OUTPUT_DIR}}/fastq > {log.mkdir} 2>&1 ; "
        #" for FILE in {input.forward_fastqs}; "
        #" do "
        #"       BASE_FILENAME=`basename ${{FILE}}`; "
        #"       ln -f ${{FILE}} ${{OUTPUT_DIR}}/fastq/${{BASE_FILENAME%{params.fastq_extensions}}}_R1_001{params.fastq_extensions} >> {log.ln} 2>&1; "
        #" done; "
        #" for FILE in {input.reverse_fastqs}; "
        #" do "
        #"       BASE_FILENAME=`basename ${{FILE}}`; "
        #"       ln -f ${{FILE}} ${{OUTPUT_DIR}}/fastq/${{BASE_FILENAME%{params.fastq_extensions}}}_R2_001{params.fastq_extensions} >> {log.ln} 2>&1; "
        #" done; "

rule threeddna: #
    input:
        fasta=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..threeddna_" + wildcards.hic_scaffolding_parameters]["prev_stage"],
                                                                                  wildcards.prev_stage_parameters,
                                                                                  wildcards.genome_prefix,
                                                                                  wildcards.haplotype),
        merged_nodups=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..threeddna_{hic_scaffolding_parameters}/{haplotype}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.merged_nodups.txt"
    params:
        restriction_seq=config["hic_enzyme_set"]  if config["hic_enzyme_set"] not in config["no_motif_enzyme_sets"] else "none",
        fastq_extensions=config["fastq_extension"]
    output:
        draft_fasta=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..threeddna_{hic_scaffolding_parameters, [^/]+}/{haplotype, [^.]+}/scaffolding/{genome_prefix, [^/]+}.input.{haplotype}.fasta",
        rawchrom_hic=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..threeddna_{hic_scaffolding_parameters, [^/]+}/{haplotype, [^.]+}/scaffolding/{genome_prefix, [^/]+}.input.{haplotype}.rawchrom.hic",
        rawchrom_assembly=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..threeddna_{hic_scaffolding_parameters, [^/]+}/{haplotype, [^.]+}/scaffolding/{genome_prefix, [^/]+}.input.{haplotype}.rawchrom.assembly",
        hic_fasta=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..threeddna_{hic_scaffolding_parameters, [^/]+}/{haplotype, [^.]+}/scaffolding/{genome_prefix, [^/]+}.input.{haplotype}_HiC.fasta",
        alias_fasta=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..threeddna_{hic_scaffolding_parameters, [^/]+}/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype, [^.]+}.fasta",

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
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("threeddna"),
        cpus=parameters["threads"]["threeddna"] ,
        time=parameters["time"]["threeddna"],
        mem=parameters["memory_mb"]["threeddna"]
    threads: parameters["threads"]["threeddna"]

    shell:
        " OUTPUT_DIR=`dirname {output.draft_fasta}`; "
        " SCRIPT=`realpath -s workflow/external_tools/3d-dna/run-3ddna-pipeline.sh`; "
        " INPUT_FASTA=`realpath -s {input.fasta}`; "
        " THREEDDNA_LOG=`realpath -s  {log.threeddna}`; "
        " LN_LOG=`realpath -s {log.ln}`; "
        " > ${{LN_LOG}}; "
        " ln -s ${{INPUT_FASTA}} {output.draft_fasta} >> ${{LN_LOG}} 2>&1; "
        " cd ${{OUTPUT_DIR}}; "
        " ${{SCRIPT}} `basename {output.draft_fasta}` `basename {input.merged_nodups}` > ${{THREEDDNA_LOG}} 2>&1; "
        " ln -s {wildcards.haplotype}/scaffolding/{wildcards.genome_prefix}.input.{wildcards.haplotype}_HiC.fasta "
        " ../../`basename {output.alias_fasta}` >> ${{LN_LOG}} 2>&1"
