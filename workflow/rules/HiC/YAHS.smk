localrules: create_links_for_yahs_files
#ruleorder: create_links_for_yahs_files > juicer_tools_pre
rule yahs: #
    input:
        #hic_qc=lambda wildcards: expand(out_dir_path / "{0}/{1}/{2}/alignment/{3}/{4}.{0}.{3}.{2}.rmdup.inter_{5}.txt".format(stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..yahs_" + wildcards.hic_scaffolding_parameters]["prev_stage"],
        #                                                                                                                      wildcards.prev_stage_parameters,
        #                                                                                                                      wildcards.haplotype,
        #                                                                                                                      stage_dict["hic_scaffolding"]["parameters"][wildcards.prev_stage_parameters + "..yahs_" + wildcards.hic_scaffolding_parameters]["option_set"]["phasing_kmer_length"],
        #                                                                                                                      wildcards.genome_prefix, "{mapq}"),
        #                                mapq=parameters["tool_options"]["juicer_tools_qc"]["mapping_quality"],
        #                                allow_missing=True),
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
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("yahs"),
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
        links_bed=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..yahs_{hic_scaffolding_parameters, [^/]+}/{haplotype, [^./]+}/scaffolding/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype}.bed",
        liftover_agp=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..yahs_{hic_scaffolding_parameters, [^/]+}/{haplotype, [^./]+}/scaffolding/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype}.liftover.agp",
        assembly=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..yahs_{hic_scaffolding_parameters, [^/]+}/{haplotype, [^./]+}/scaffolding/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype}.assembly",
        assembly_agp=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..yahs_{hic_scaffolding_parameters, [^/]+}/{haplotype, [^.]+}/scaffolding/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype}.assembly.agp",
        log=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..yahs_{hic_scaffolding_parameters, [^/]+}/{haplotype, [^./]+}/scaffolding/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype}.log",
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
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("yahs_juicer_pre"),
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
        hic=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype, [^.]+}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.hic",
    params:
        resolution_list=" ".join(map(str,parameters["tool_options"]["juicer_tools_pre"]["resolution_list"]))
    log:
        grep0=output_dict["log"]  / "juicer_tools_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.grep0.log",
        sed0=output_dict["log"]  / "juicer_tools_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.sed0.log",
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
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("juicer_tools_pre"),
        cpus=parameters["threads"]["juicer_tools_pre"] ,
        time=parameters["time"]["juicer_tools_pre"],
        mem=parameters["memory_mb"]["juicer_tools_pre"]
    threads: parameters["threads"]["juicer_tools_pre"]

    shell: # juicer_tools elder than 1.9.9 seems to be incompartible with yahs
        " RESOLUTION_LIST=({params.resolution_list}); "
        " RESOLUTION_LIST_LEN=${{#RESOLUTION_LIST[@]}}; "
        " SCALE=`grep 'scale factor:' {input.yahs_juicer_pre_log} 2>{log.grep0} | sed 's/.*scale factor: //' 2>{log.sed0}`; "
        " RESOLUTION_OPTION_LIST=$(( RESOLUTION_LIST[0]/SCALE )); "
        " for (( i=1; i<$RESOLUTION_LIST_LEN; i++ )); "
        "   do"
        "   RESOLUTION_OPTION_LIST=$RESOLUTION_OPTION_LIST\",\"$(( RESOLUTION_LIST[$i]/SCALE ));"
        "   done;"
        " java -jar -Xmx{resources.mem}m workflow/external_tools/juicer/juicer_tools.1.9.9_jcuda.0.8.jar pre "
        " -r ${{RESOLUTION_OPTION_LIST}} " #--threads {threads}  
        " {input.yahs_juicer_pre_bed} {output.hic} <(cat {input.yahs_juicer_pre_log} 2>{log.cat} | "
        " grep PRE_C_SIZE 2>{log.grep} | awk '{{print $2\" \"$3}}' 2>{log.awk}) > {log.juicer} 2>&1; "
        #
        #" java -jar -Xmx{resources.mem}m workflow/external_tools/juicer/juicer_tools.1.9.9_jcuda.0.8.jar pre " #--threads {threads}
        #" {input.yahs_juicer_pre_bed} {output.hic} <(cat {input.yahs_juicer_pre_log} 2>{log.cat} | "
        #" grep PRE_C_SIZE 2>{log.grep} | awk '{{print $2\" \"$3}}' 2>{log.awk}) > {log.juicer} 2>&1; "


rule create_links_for_yahs_files: #
    input:
        hic=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.hic",
        liftover_agp=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.liftover.agp",
        assembly=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.assembly",
        assembly_agp=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{haplotype}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.assembly.agp",
    output:
        hic=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..yahs_{hic_scaffolding_parameters, [^/]+}/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype, [^./]+}.hic",
        liftover_agp=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..yahs_{hic_scaffolding_parameters, [^/]+}/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype, [^./]+}.liftover.agp",
        assembly=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..yahs_{hic_scaffolding_parameters, [^/]+}/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype, [^./]+}.assembly",
        assembly_agp=out_dir_path / "hic_scaffolding/{prev_stage_parameters, [^/]+}..yahs_{hic_scaffolding_parameters, [^/]+}/{genome_prefix, [^/]+}.hic_scaffolding.{haplotype, [^./]+}.assembly.agp",
    log:
        ln=output_dict["log"]  / "create_links_for_yahs_files.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "create_links_for_yahs_files.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_links_for_yahs_files.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_links_for_yahs_files.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("create_links_for_yahs_files"),
        cpus=parameters["threads"]["create_fastq_links"] ,
        time=parameters["time"]["create_fastq_links"],
        mem=parameters["memory_mb"]["create_fastq_links"]
    threads: parameters["threads"]["create_fastq_links"]

    shell:
        " ln -sf {wildcards.haplotype}/scaffolding/`basename {input.hic}` {output.hic} > {log.ln} 2>&1; "
        " ln -sf {wildcards.haplotype}/scaffolding/`basename {input.liftover_agp}` {output.liftover_agp} >> {log.ln} 2>&1; "
        " ln -sf {wildcards.haplotype}/scaffolding/`basename {input.assembly}` {output.assembly} >> {log.ln} 2>&1; "
        " ln -sf {wildcards.haplotype}/scaffolding/`basename {input.assembly_agp}` {output.assembly_agp} >> {log.ln} 2>&1; "

rule extract_yahs_contigs: #
    input:
        liftover_agp=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{genome_prefix}.hic_scaffolding.{haplotype}.liftover.agp",
        original_contigs=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta" % (stage_dict["hic_scaffolding"]["prev_stage"],
                                                                                                              stage_dict["hic_scaffolding"]["prev_stage"])) ,
        original_contigs_fai=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta.fai" % (stage_dict["hic_scaffolding"]["prev_stage"],
                                                                                                                      stage_dict["hic_scaffolding"]["prev_stage"]))
    output:
        new_contigs_fasta=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{genome_prefix}.hic_scaffolding.{haplotype}.new_contigs.fasta",
    log:
        std=output_dict["log"]  / "extract_yahs_contigs.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.std.log",
        cluster_log=output_dict["cluster_log"] / "extract_yahs_contigs.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "extract_yahs_contigs.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "extract_yahs_contigs.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("extract_yahs_contigs"),
        cpus=parameters["threads"]["extract_yahs_contigs"] ,
        time=parameters["time"]["extract_yahs_contigs"],
        mem=parameters["memory_mb"]["extract_yahs_contigs"]
    threads: parameters["threads"]["extract_yahs_contigs"]

    shell:
        " agp_to_fasta -o {output.new_contigs_fasta} {input.liftover_agp} {input.original_contigs} > {log.std} 2>&1; "

rule create_transfer_agp: #
    input:
        liftover_agp=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{genome_prefix}.hic_scaffolding.{haplotype}.liftover.agp",
        assembly_agp=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{genome_prefix}.hic_scaffolding.{haplotype}.assembly.agp",
    output:
        transfer_agp=out_dir_path / "hic_scaffolding/{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}/{genome_prefix}.hic_scaffolding.{haplotype}.transfer.agp"
    log:
        std=output_dict["log"]  / "create_transfer_agp.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.std.log",
        cluster_log=output_dict["cluster_log"] / "create_transfer_agp.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_transfer_agp.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_transfer_agp.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("create_transfer_agp"),
        cpus=parameters["threads"]["extract_yahs_contigs"] ,
        time=parameters["time"]["extract_yahs_contigs"],
        mem=parameters["memory_mb"]["extract_yahs_contigs"]
    threads: parameters["threads"]["extract_yahs_contigs"]

    shell:
        " ./workflow/scripts/hic_scaffolding/create_transfer_agp.py "
        " -a {input.assembly_agp} -l {input.liftover_agp} -o {output.transfer_agp} > {log.std} 2>&1; "