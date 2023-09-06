localrules: generate_prescaffolding_agp


rule generate_site_positions: #
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta"
    output:
        restriction_site_file=out_dir_path / ("{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}_%s.txt" % config["hic_enzyme_set"])
    params:
        restriction_seq=config["hic_enzyme_set"]
    log:
        sites=output_dict["log"]  / "generate_site_positions.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.sites.log",
        cluster_log=output_dict["cluster_log"] / "generate_site_positions.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "generate_site_positions.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "generate_site_positions.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("generate_site_positions"),
        cpus=parameters["threads"]["generate_site_positions"] ,
        time=parameters["time"]["generate_site_positions"],
        mem=parameters["memory_mb"]["generate_site_positions"]
    threads: parameters["threads"]["generate_site_positions"]

    shell:
        #" ln -sf ../../../../../../{input.fasta} {output.alias_fasta} > {log.ln} 2>&1; "
        " OUTPUT_PREFIX={output.restriction_site_file}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%_{params.restriction_seq}.txt}}; "
        " ./workflow/external_tools/juicer/misc/generate_site_positions.py {params.restriction_seq} ${{OUTPUT_PREFIX}} {input.fasta} > {log.sites} 2>&1; "

rule generate_prescaffolding_agp:
    input:
        reference_fai="{prefix}.fasta.fai",
    output:
        pre_agp="{prefix}.pre.agp"
    log:
        std="{prefix}.generatepre_scaffolding_agp.log",
        cluster_log="{prefix}.generate_prescaffolding_agp.cluster.log",
        cluster_err="{prefix}.generate_prescaffolding_agp.cluster.err"
    benchmark:
        "{prefix}.generate_prescaffolding_agp.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("generate_prescaffolding_agp"),
        cpus=parameters["threads"]["generate_prescaffolding_agp"] ,
        time=parameters["time"]["generate_prescaffolding_agp"],
        mem=parameters["memory_mb"]["generate_prescaffolding_agp"]
    threads: parameters["threads"]["generate_prescaffolding_agp"]

    shell:
        " ./workflow/scripts/hic_scaffolding/generate_agp_from_fai.py -f {input.reference_fai} -o {output.pre_agp} > {log.std} 2>&1;"

rule yahs_juicer_pre_prescaffolding: #
    input:
        reference_fai=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta.fai",
        bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
        bai=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam.bai",
        agp=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.pre.agp"
    output:
        links_bed=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.pre.bed",
        liftover_agp=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.pre.liftover.agp",
        liftover_syn=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.pre.liftover.syn",
        assembly=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.pre.assembly",
        assembly_original=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.pre.original.assembly",
        assembly_agp=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.pre.assembly.agp",
        log=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.pre.yahs.juicer_pre.log",
    log:
        #std=output_dict["log"]  / "yahs_juicer_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.log",
        mv=output_dict["log"]  / "yahs_juicer_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.mv.log",
        cut=output_dict["log"]  / "yahs_juicer_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.cut.log",
        rename=output_dict["log"]  / "yahs_juicer_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.rename.log",
        cluster_log=output_dict["cluster_log"] / "yahs_juicer_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "yahs_juicer_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "yahs_juicer_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("yahs_juicer_pre_prescaffolding"),
        cpus=parameters["threads"]["yahs_juicer_pre"] ,
        time=parameters["time"]["yahs_juicer_pre"],
        mem=parameters["memory_mb"]["yahs_juicer_pre"]
    threads: parameters["threads"]["yahs_juicer_pre"]
    shell:
        " OUTPUT_PREFIX={output.links_bed}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.bed}}; "
        " juicer pre -a -o ${{OUTPUT_PREFIX}} {input.bam} {input.agp} {input.reference_fai} > {output.log} 2>&1;"
        " mv ${{OUTPUT_PREFIX}}.txt {output.links_bed} > {log.mv} 2>&1; "
        " cut -f 1,6 {output.liftover_agp} > {output.liftover_syn} 2>{log.cut}; "
        " ./workflow/scripts/hic_scaffolding/rename_scaffolds_in_assembly_file.py -a {output.assembly} "
        " --scaffold_syn_file {output.liftover_syn} --syn_file_key_column 0 --syn_file_value_column 1 "
        " -o {output.assembly_original} > {log.rename} 2>&1; "

rule juicer_tools_pre_prescaffolding: #
    input:
        yahs_juicer_pre_log=rules.yahs_juicer_pre_prescaffolding.output.log,
        yahs_juicer_pre_bed=rules.yahs_juicer_pre_prescaffolding.output.links_bed
    output:
        hic=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.pre.hic",
    params:
        resolution_list=" ".join(map(str,parameters["tool_options"]["juicer_tools_pre"]["resolution_list"]))
    log:
        grep0=output_dict["log"]  / "juicer_tools_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.grep0.log",
        sed0=output_dict["log"]  / "juicer_tools_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.sed0.log",
        juicer=output_dict["log"]  / "juicer_tools_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.juicer.log",
        cat=output_dict["log"]  / "juicer_tools_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.cat.log",
        grep=output_dict["log"]  / "juicer_tools_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.grep.log",
        awk=output_dict["log"]  / "juicer_tools_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.awk.log",
        cluster_log=output_dict["cluster_log"] / "juicer_tools_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "juicer_tools_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "juicer_tools_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
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
        "   do"
        "   RESOLUTION_OPTION_LIST=$RESOLUTION_OPTION_LIST\",\"$(( RESOLUTION_LIST[$i]/SCALE ));"
        "   done;"
        " java -jar -Xmx{resources.mem}m workflow/external_tools/juicer/juicer_tools.1.9.9_jcuda.0.8.jar pre "
        " -r ${{RESOLUTION_OPTION_LIST}} " #--threads {threads}  
        " {input.yahs_juicer_pre_bed} {output.hic} <(cat {input.yahs_juicer_pre_log} 2>{log.cat} | "
        " grep PRE_C_SIZE 2>{log.grep} | awk '{{print $2\" \"$3}}' 2>{log.awk}) > {log.juicer} 2>&1; "

rule bam2bed:
    input:
        bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam"
    output:
        #bed=out_dir_path  / ("{assembly_stage}/{assembler}/{haplotype}/alignment/%s.{assembly_stage}.{assembler}.{haplotype}.bwa.filtered.rmdup.bed"  % config["genome_name"]),
        bed=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bed"
    log:
        convert=output_dict["log"] / "bam2bed.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.convert.log",
        sort=output_dict["log"] / "bam2bed.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.sort.log",
        cluster_log=output_dict["cluster_log"] / "bam2bed.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bam2bed.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bam2bed.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("bam_to_bed"),
        cpus=parameters["threads"]["bam2bed"] ,
        time=parameters["time"]["bam2bed"],
        mem=parameters["memory_mb"]["bam2bed"]
    threads: parameters["threads"]["bam2bed"]

    shell:
        " bamToBed -i {input.bam} 2>{log.convert} | sort -S{resources.mem}M -k 4 > {output.bed} 2>{log.sort}"
