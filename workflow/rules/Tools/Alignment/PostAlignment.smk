localrules: generate_prescaffolding_agp

rule generate_site_positions: #
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/")
    output:
        restriction_site_file="{fasta_dir}/{fasta_prefix}_%s.txt" % config["hic_enzyme_set"]
    params:
        restriction_seq=config["hic_enzyme_set"]
    log:
        sites="{fasta_dir}/log/generate_site_positions.{fasta_prefix}.sites.log",
        cluster_log="{fasta_dir}/log/generate_site_positions.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/log/generate_site_positions.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/log/generate_site_positions.{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("generate_site_positions"),
        cpus=parameters["threads"]["generate_site_positions"] ,
        time=parameters["time"]["generate_site_positions"],
        mem=parameters["memory_mb"]["generate_site_positions"]
    threads: parameters["threads"]["generate_site_positions"]

    shell:
        " OUTPUT_PREFIX={output.restriction_site_file}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%_{params.restriction_seq}.txt}}; "
        " ./workflow/external_tools/juicer/misc/generate_site_positions.py {params.restriction_seq} ${{OUTPUT_PREFIX}} {input.fasta} > {log.sites} 2>&1; "

rule generate_prescaffolding_agp:
    input:
        reference_fai="{fasta_dir}/{fasta_prefix}.fasta.fai",
        log_dir=ancient("{fasta_dir}/log/")
    output:
        pre_agp="{fasta_dir}/{fasta_prefix}.pre.agp"
    log:
        std="{fasta_dir}/log/generatepre_scaffolding_agp.{fasta_prefix}.log",
        cluster_log="{fasta_dir}/log/generate_prescaffolding_agp.cluster.{fasta_prefix}.log",
        cluster_err="{fasta_dir}/log/generate_prescaffolding_agp.cluster.{fasta_prefix}.err"
    benchmark:
        "{fasta_dir}/log/generate_prescaffolding_agp.benchmark.{fasta_prefix}.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("generate_prescaffolding_agp"),
        cpus=parameters["threads"]["generate_prescaffolding_agp"] ,
        time=parameters["time"]["generate_prescaffolding_agp"],
        mem=parameters["memory_mb"]["generate_prescaffolding_agp"]
    threads: parameters["threads"]["generate_prescaffolding_agp"]

    shell:
        " ./workflow/scripts/hic_scaffolding/generate_agp_from_fai.py -f {input.reference_fai} -o {output.pre_agp} > {log.std} 2>&1;"


rule yahs_juicer_pre_prescaffolding: #
    input:
        reference_fai="{fasta_dir}/{fasta_prefix}.fasta.fai",
        bam="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.bam",
        csi="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.bam.csi",
        agp="{fasta_dir}/{fasta_prefix}.pre.agp",
        log_dir=ancient("{fasta_dir}/log/")
    output:
        links_bed="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.pre.mapq{min_mapq}.bed",
        liftover_agp="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.pre.mapq{min_mapq}.liftover.agp",
        liftover_syn="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.pre.mapq{min_mapq}.liftover.syn",
        assembly="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.pre.mapq{min_mapq}.assembly",
        assembly_original="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.pre.mapq{min_mapq}.original.assembly",
        assembly_agp="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.pre.mapq{min_mapq}.assembly.agp",
        log="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.pre.mapq{min_mapq}.yahs.juicer_pre.log",

    log:
        mv="{fasta_dir}/log/yahs_juicer_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.mv.log",
        cut="{fasta_dir}/log/yahs_juicer_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.cut.log",
        rename="{fasta_dir}/log/yahs_juicer_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.rename.log",
        cluster_log="{fasta_dir}/log/yahs_juicer_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.cluster.log",
        cluster_err="{fasta_dir}/log/yahs_juicer_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.cluster.err"
    benchmark:
        "{fasta_dir}/log/yahs_juicer_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("yahs_juicer_pre_prescaffolding"),
        cpus=parameters["threads"]["yahs_juicer_pre"] ,
        time=parameters["time"]["yahs_juicer_pre"],
        mem=parameters["memory_mb"]["yahs_juicer_pre"]
    threads: parameters["threads"]["yahs_juicer_pre"]
    shell:
        " OUTPUT_PREFIX={output.links_bed}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.bed}}; "
        " juicer pre -a -q {wildcards.min_mapq} -o ${{OUTPUT_PREFIX}} {input.bam} {input.agp} {input.reference_fai} > {output.log} 2>&1;"
        " mv ${{OUTPUT_PREFIX}}.txt {output.links_bed} > {log.mv} 2>&1; "
        " cut -f 1,6 {output.liftover_agp} > {output.liftover_syn} 2>{log.cut}; "
        " ./workflow/scripts/hic_scaffolding/rename_scaffolds_in_assembly_file.py -a {output.assembly} "
        "     --scaffold_syn_file {output.liftover_syn} --syn_file_key_column 0 --syn_file_value_column 1 "
        "     -o {output.assembly_original} > {log.rename} 2>&1; "

rule juicer_tools_pre_prescaffolding: #
    input:
        yahs_juicer_pre_log=rules.yahs_juicer_pre_prescaffolding.output.log,
        yahs_juicer_pre_bed=rules.yahs_juicer_pre_prescaffolding.output.links_bed,
        log_dir=ancient("{fasta_dir}/log/")
    output:
        hic="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.rmdup.pre.mapq{min_mapq}.hic",
    params:
        resolution_list=" ".join(map(str, parameters["tool_options"]["juicer_tools_pre"]["resolution_list"]))
    log:
        grep0="{fasta_dir}/log/juicer_tools_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.grep0.log",
        sed0="{fasta_dir}/log/juicer_tools_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.sed0.log",
        juicer="{fasta_dir}/log/juicer_tools_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.juicer.log",
        cat="{fasta_dir}/log/juicer_tools_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.cat.log",
        grep="{fasta_dir}/log/juicer_tools_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.grep.log",
        awk="{fasta_dir}/log/juicer_tools_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.awk.log",
        cluster_log="{fasta_dir}/log/juicer_tools_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.cluster.log",
        cluster_err="{fasta_dir}/log/juicer_tools_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.cluster.err"
    benchmark:
        "{fasta_dir}/log/juicer_tools_pre_prescaffolding.{fasta_prefix}.{phasing_kmer_length}.mapq{min_mapq}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("juicer_tools_pre_prescaffolding"),
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
        " do"
        "     RESOLUTION_OPTION_LIST=$RESOLUTION_OPTION_LIST\",\"$(( RESOLUTION_LIST[$i]/SCALE ));"
        " done;"
        " java -jar -Xmx{resources.mem}m workflow/external_tools/juicer/juicer_tools.1.9.9_jcuda.0.8.jar pre "
        "     -r ${{RESOLUTION_OPTION_LIST}} {input.yahs_juicer_pre_bed} {output.hic} <(cat {input.yahs_juicer_pre_log} 2>{log.cat} | "
        "     grep PRE_C_SIZE 2>{log.grep} | awk '{{print $2\" \"$3}}' 2>{log.awk}) > {log.juicer} 2>&1; "

