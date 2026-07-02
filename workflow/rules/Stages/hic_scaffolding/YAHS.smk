#localrules: create_links_for_yahs_files
#ruleorder: yahs_juicer_pre > get_seq_len

rule yahs: #
    input:
        bam=lambda wildcards: config["out_dir"]  / "{0}/{1}/{4}.{0}.{2}/alignment/{3}/{4}.{0}.{2}.{3}.rmdup.bam".format(stage_dict["hic_scaffolding"].prev_stage,
                                                                                                                        get_prev_stage_parameters(wildcards.parameters),
                                                                                                                        wildcards.haplotype,
                                                                                                                        stage_dict["hic_scaffolding"].parameters[wildcards.parameters]["option_set"]["phasing_kmer_length"],
                                                                                                                        wildcards.genome_prefix),
        bam_stats=lambda wildcards: config["out_dir"]  / "{0}/{1}/{4}.{0}.{2}/alignment/{3}/{4}.{0}.{2}.{3}.rmdup.bam.general_stats".format(stage_dict["hic_scaffolding"].prev_stage,
                                                                                                                                            get_prev_stage_parameters(wildcards.parameters),
                                                                                                                                            wildcards.haplotype,
                                                                                                                                            stage_dict["hic_scaffolding"].parameters[wildcards.parameters]["option_set"]["phasing_kmer_length"],
                                                                                                                                            wildcards.genome_prefix),
        reference=lambda wildcards: config["out_dir"]  / ("%s/%s/%s.%s.%s.fasta" % (stage_dict["hic_scaffolding"].prev_stage,
                                                                                    get_prev_stage_parameters(wildcards.parameters),
                                                                                    wildcards.genome_prefix,
                                                                                    stage_dict["hic_scaffolding"].prev_stage,
                                                                                    wildcards.haplotype)),
        reference_fai=lambda wildcards: config["out_dir"]  / ("%s/%s/%s.%s.%s.fasta.fai" % (stage_dict["hic_scaffolding"].prev_stage,
                                                                                            get_prev_stage_parameters(wildcards.parameters),
                                                                                            wildcards.genome_prefix,
                                                                                            stage_dict["hic_scaffolding"].prev_stage,
                                                                                            wildcards.haplotype)),
    output:
        fasta=config["out_dir"]  / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^.]+}/scaffolding/{genome_prefix}_scaffolds_final.fa",
        bin=config["out_dir"]  / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^.]+}/scaffolding/{genome_prefix}.bin",
        agp=config["out_dir"] / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^.]+}/scaffolding/{genome_prefix}_scaffolds_final.agp",
        alias=config["out_dir"] / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^/]+}.fasta",
    params:
        min_contig_len=lambda wildcards: stage_dict["hic_scaffolding"].parameters[wildcards.parameters]["option_set"]["min_contig_len"],
        min_mapping_quality=lambda wildcards: stage_dict["hic_scaffolding"].parameters[wildcards.parameters]["option_set"]["min_mapping_quality"],
        restriction_seq=parse_option(config["hic_enzyme_set"], parameters["tool_options"]["yahs"]["restriction_seq"], "-e",)
    log:
        yahs=config["out_dir"]  / "log/yahs.{parameters}.{genome_prefix}.hic_scaffolding.{haplotype}.yahs.log",
        ln=config["out_dir"] / "log/yahs.{parameters}.{genome_prefix}.hic_scaffolding.{haplotype}.ln.log",
        cluster_log=config["out_dir"] / "log/yahs.{parameters}.{genome_prefix}.hic_scaffolding.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/yahs.{parameters}.{genome_prefix}.hic_scaffolding.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"]  / "log/yahs.{parameters}.{genome_prefix}.hic_scaffolding.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("yahs"),
        cpus=parameters["threads"]["yahs"] ,
        time=parameters["time"]["yahs"],
        mem=parameters["memory_mb"]["yahs"]
    threads: parameters["threads"]["yahs"]

    shell:
        " OUTPUT_PREFIX={output.fasta}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%_scaffolds_final.fa}}; "
        " yahs -q {params.min_mapping_quality} -l {params.min_contig_len} -o ${{OUTPUT_PREFIX}} "
        "     {input.reference} {input.bam} > {log.yahs} 2>&1;"
        " ln -sf {wildcards.genome_prefix}.hic_scaffolding.{wildcards.haplotype}/scaffolding/`basename {output.fasta}` {output.alias} > {log.ln} 2>&1"


rule yahs_juicer_pre: #
    input:
        reference_fai=rules.yahs.input.reference_fai,
        bin=rules.yahs.output.bin,
        agp=rules.yahs.output.agp,
        log_dir=ancient(config["out_dir"]  / "hic_scaffolding/{parameters}/log/")
    output:
        links_bed=config["out_dir"]  / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^.]+}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.bed",
        liftover_agp=config["out_dir"]  / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^.]+}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.liftover.agp",
        assembly=config["out_dir"]  / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^.]+}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.assembly",
        assembly_agp=config["out_dir"]  / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^.]+}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.assembly.agp",
        log=config["out_dir"]  / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^.]+}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.log",
    log:
        mv=config["out_dir"]  / "hic_scaffolding/{parameters}/log/yahs_juicer_pre.{parameters}.{genome_prefix}.hic.scaffolding.{haplotype}.mv.log",
        cluster_log=config["out_dir"]  / "hic_scaffolding/{parameters}/log/yahs_juicer_pre.{parameters}.{genome_prefix}.hic.scaffolding.{haplotype}.cluster.log",
        cluster_err=config["out_dir"]  / "hic_scaffolding/{parameters}/log/yahs_juicer_pre.{parameters}.{genome_prefix}.hic.scaffolding.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"]  / "hic_scaffolding/{parameters}/log/yahs_juicer_pre.{genome_prefix}.hic.scaffolding.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("yahs_juicer_pre"),
        cpus=parameters["threads"]["yahs_juicer_pre"] ,
        time=parameters["time"]["yahs_juicer_pre"],
        mem=parameters["memory_mb"]["yahs_juicer_pre"]
    threads: parameters["threads"]["yahs_juicer_pre"]

    shell:
        " OUTPUT_PREFIX={output.links_bed}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.bed}}; "
        " juicer pre -a -o ${{OUTPUT_PREFIX}} {input.bin} {input.agp} {input.reference_fai} > {output.log} 2>&1; "
        " mv ${{OUTPUT_PREFIX}}.txt {output.links_bed} > {log.mv} 2>&1; "

rule juicer_tools_pre: #
    input:
        yahs_juicer_pre_log=rules.yahs_juicer_pre.output.log,
        yahs_juicer_pre_bed=rules.yahs_juicer_pre.output.links_bed,
        log_dir=ancient(config["out_dir"]  / "hic_scaffolding/{parameters}/log/")
    output:
        hic=config["out_dir"]  / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, hap[^/]+}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.hic",

    params:
        resolution_list=" ".join(map(str, parameters["tool_options"]["juicer_tools_pre"]["resolution_list"]))
    log:
        grep0=config["out_dir"]  / "hic_scaffolding/{parameters}/log/yahs_juicer_pre.{parameters}.{genome_prefix}.hic.scaffolding.{haplotype}.grep0.log",
        sed0=config["out_dir"]  / "hic_scaffolding/{parameters}/log/yahs_juicer_pre.{parameters}.{genome_prefix}.hic.scaffolding.{haplotype}.sed0.log",
        juicer=config["out_dir"]  / "hic_scaffolding/{parameters}/log/yahs_juicer_pre.{parameters}.{genome_prefix}.hic.scaffolding.{haplotype}.juicer.log",
        cat=config["out_dir"]  / "hic_scaffolding/{parameters}/log/yahs_juicer_pre.{parameters}.{genome_prefix}.hic.scaffolding.{haplotype}.cat.log",
        grep=config["out_dir"]  / "hic_scaffolding/{parameters}/log/yahs_juicer_pre.{parameters}.{genome_prefix}.hic.scaffolding.{haplotype}.grep.log",
        awk=config["out_dir"]  / "hic_scaffolding/{parameters}/log/yahs_juicer_pre.{parameters}.{genome_prefix}.hic.scaffolding.{haplotype}.awk.log",
        cluster_log=config["out_dir"]  / "hic_scaffolding/{parameters}/log/yahs_juicer_pre.{parameters}.{genome_prefix}.hic.scaffolding.{haplotype}.cluster.log",
        cluster_err=config["out_dir"]  / "hic_scaffolding/{parameters}/log/yahs_juicer_pre.{parameters}.{genome_prefix}.hic.scaffolding.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"]  / "hic_scaffolding/{parameters}/log/yahs_juicer_pre.{parameters}.{genome_prefix}.hic.scaffolding.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
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
        " do "
        "     RESOLUTION_OPTION_LIST=$RESOLUTION_OPTION_LIST\",\"$(( RESOLUTION_LIST[$i]/SCALE )); "
        " done;"
        " java -jar -Xmx{resources.mem}m workflow/external_tools/juicer/juicer_tools.1.9.9_jcuda.0.8.jar pre "
        "     -r ${{RESOLUTION_OPTION_LIST}} {input.yahs_juicer_pre_bed} {output.hic} <(cat {input.yahs_juicer_pre_log} 2>{log.cat} | "
        "     grep PRE_C_SIZE 2>{log.grep} | awk '{{print $2\" \"$3}}' 2>{log.awk}) > {log.juicer} 2>&1; "

use rule create_local_links as create_links_for_yahs_files with: #
    input:
        hic=config["out_dir"]  / "hic_scaffolding/{parameters}/{genome_prefix}.hic_scaffolding.{haplotype}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.hic",
        liftover_agp=config["out_dir"]  / "hic_scaffolding/{parameters}/{genome_prefix}.hic_scaffolding.{haplotype}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.liftover.agp",
        assembly=config["out_dir"]  / "hic_scaffolding/{parameters}/{genome_prefix}.hic_scaffolding.{haplotype}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.assembly",
        assembly_agp=config["out_dir"]  / "hic_scaffolding/{parameters}/{genome_prefix}.hic_scaffolding.{haplotype}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.assembly.agp",
        log_dir=ancient(config["out_dir"]  / "hic_scaffolding/{parameters}/log/")
    output:
        hic=config["out_dir"]  / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype}.hic",
        liftover_agp=config["out_dir"]  / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype}.liftover.agp",
        assembly=config["out_dir"]  / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype}.assembly",
        assembly_agp=config["out_dir"]  / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype}.assembly.agp",
    log:
        ln=config["out_dir"]  / "hic_scaffolding/{parameters}/log/create_links_for_yahs_files.{parameters}.{genome_prefix}.hic_scaffolding.{haplotype}.ln.log",


rule extract_yahs_contigs: #
    input:
        liftover_agp=config["out_dir"] / "hic_scaffolding/{parameters}/{genome_prefix}.hic_scaffolding.{haplotype}.liftover.agp",
        original_contigs=lambda wildcards: config["out_dir"]  / ("%s/%s/%s.%s.%s.fasta" % (stage_dict["hic_scaffolding"].prev_stage,
                                                                                           get_prev_stage_parameters(wildcards.parameters),
                                                                                           wildcards.genome_prefix,
                                                                                           stage_dict["hic_scaffolding"].prev_stage,
                                                                                           wildcards.haplotype)),
        original_contigs_fai=lambda wildcards: config["out_dir"]  / ("%s/%s/%s.%s.%s.fasta.fai" % (stage_dict["hic_scaffolding"].prev_stage,
                                                                                                   get_prev_stage_parameters(wildcards.parameters),
                                                                                                   wildcards.genome_prefix,
                                                                                                   stage_dict["hic_scaffolding"].prev_stage,
                                                                                                   wildcards.haplotype)),
        log_dir=ancient(config["out_dir"]  / "hic_scaffolding/{parameters}/log/")
    output:
        new_contigs_fasta=config["out_dir"] / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, [^/]+}.new_contigs.fasta",
    log:
        std=config["out_dir"]  / "hic_scaffolding/{parameters}/log/extract_yahs_contigs.{parameters}.{genome_prefix}.{haplotype}.std.log",
        cluster_log=config["out_dir"]  / "hic_scaffolding/{parameters}/log/extract_yahs_contigs.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"]  / "hic_scaffolding/{parameters}/log/extract_yahs_contigs.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"]  / "hic_scaffolding/{parameters}/log/extract_yahs_contigs.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("extract_yahs_contigs"),
        cpus=parameters["threads"]["extract_yahs_contigs"] ,
        time=parameters["time"]["extract_yahs_contigs"],
        mem=parameters["memory_mb"]["extract_yahs_contigs"]
    threads: parameters["threads"]["extract_yahs_contigs"]

    shell:
        " agp_to_fasta -o {output.new_contigs_fasta} {input.liftover_agp} {input.original_contigs} > {log.std} 2>&1; "

rule create_transfer_agp: #
    input:
        liftover_agp=config["out_dir"] / "hic_scaffolding/{parameters}/{genome_prefix}.hic_scaffolding.{haplotype}.liftover.agp",
        assembly_agp=config["out_dir"] / "hic_scaffolding/{parameters}/{genome_prefix}.hic_scaffolding.{haplotype}.assembly.agp",
        log_dir=ancient(config["out_dir"]  / "hic_scaffolding/{parameters}/log/"),
    output:
        transfer_agp=config["out_dir"] / "hic_scaffolding/{parameters, [^/]*yahs[^/]*}/{genome_prefix}.hic_scaffolding.{haplotype, [^/]+}.transfer.agp"
    log:
        std=config["out_dir"]  / "hic_scaffolding/{parameters}/log/create_transfer_agp.{parameters}.{genome_prefix}.{haplotype}.std.log",
        cluster_log=config["out_dir"]  / "hic_scaffolding/{parameters}/log/create_transfer_agp.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"]  / "hic_scaffolding/{parameters}/log/create_transfer_agp.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"]  / "hic_scaffolding/{parameters}/log/create_transfer_agp.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_transfer_agp"),
        cpus=parameters["threads"]["extract_yahs_contigs"] ,
        time=parameters["time"]["extract_yahs_contigs"],
        mem=parameters["memory_mb"]["extract_yahs_contigs"]
    threads: parameters["threads"]["extract_yahs_contigs"]

    shell:
        " ./workflow/scripts/hic_scaffolding/create_transfer_agp.py -a {input.assembly_agp} -l {input.liftover_agp} "
        "     -o {output.transfer_agp} > {log.std} 2>&1; "
