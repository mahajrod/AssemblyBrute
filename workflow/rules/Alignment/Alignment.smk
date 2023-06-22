localrules: generate_prescaffolding_agp

rule bwa_map: #
    input:
        index=rules.bwa_index.output.index,
        reference=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
        fastq=lambda wildcards: output_dict["data"] / "fastq/hic/raw/{0}{1}".format(wildcards.fileprefix, config["fastq_extension"]) if wildcards.phasing_kmer_length == "NA" else \
                                out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}".format(config["phasing_stage"], #wildcards.assembly_stage,
                                                                                         detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
                                                                                         wildcards.haplotype,
                                                                                         wildcards.phasing_kmer_length,
                                                                                         wildcards.fileprefix,
                                                                                         config["fastq_extension"]
                                                                                         )
    output:
        bam=out_dir_path  / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{fileprefix}.bwa.bam"
    params:
        id="{0}_hic".format(config["genome_prefix"])
    log:
        map=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.map.log",
        sort=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.sort.log",
        filter=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.filter.log",
        cluster_log=output_dict["cluster_log"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["bwa_map"] ,
        time=parameters["time"]["bwa_map"],
        mem=parameters["memory_mb"]["bwa_map"]
    threads: parameters["threads"]["bwa_map"]
    shell:
        " bwa mem -SP5M -t {threads} -R  \'@RG\\tID:{params.id}\\tPU:x\\tSM:{params.id}\\tPL:illumina\\tLB:x\' "
        " {input.reference} {input.fastq} 2>{log.map} | filter_five_end.pl 2>{log.filter} | samtools view -Sb - > {output.bam} 2>{log.sort} "

rule bam_merge_pairs:
    input:
        forward_bam=lambda wildcards: out_dir_path / ("{0}/{1}/{2}/alignment/{3}/{4}.{0}.{3}.{2}.{5}{6}.bwa.bam".format(wildcards.assembly_stage,
                                                                                                                                  wildcards.parameters,
                                                                                                                                  wildcards.haplotype,
                                                                                                                                  wildcards.phasing_kmer_length,
                                                                                                                                  wildcards.genome_prefix,
                                                                                                                                  wildcards.pairprefix,
                                                                                                                                  input_forward_suffix_dict["hic"] if wildcards.phasing_kmer_length == "NA" else "_1")),
        reverse_bam=lambda wildcards: out_dir_path / ("{0}/{1}/{2}/alignment/{3}/{4}.{0}.{3}.{2}.{5}{6}.bwa.bam".format(wildcards.assembly_stage,
                                                                                                                                  wildcards.parameters,
                                                                                                                                  wildcards.haplotype,
                                                                                                                                  wildcards.phasing_kmer_length,
                                                                                                                                  wildcards.genome_prefix,
                                                                                                                                  wildcards.pairprefix,
                                                                                                                                  input_reverse_suffix_dict["hic"] if wildcards.phasing_kmer_length == "NA" else "_2")),
        reference_fai=rules.ref_faidx.output.fai
    output:
        bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{pairprefix}.bwa.bam", # TODO: make_tem
    params:
        min_mapq=parameters["tool_options"]["two_read_bam_combiner"]["mapq"],
        sort_threads=parameters["threads"]["samtools_sort"],
        sort_memory=parameters["memory_mb"]["samtools_sort"],
        #tmp_prefix=lambda wildcards: out_dir_path  / "{0}/{1}/{2}/alignment/{3}".format(wildcards.assembly_stage,
        #                                                                                wildcards.parameters,
        #                                                                                wildcards.haplotype,
        #                                                                                wildcards.pairprefix)
    log:
        merge=output_dict["log"] / "bam_merge_pairs.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.merge.log",
        view=output_dict["log"] / "bam_merge_pairs.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.view.log",
        sort=output_dict["log"] / "bam_merge_pairs.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.sort.log",
        cluster_log=output_dict["cluster_log"] / "bam_merge_pairs.{assembly_stage}.{parameters}.{phasing_kmer_length}.{genome_prefix}.{haplotype}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bam_merge_pairs.{assembly_stage}.{parameters}.{phasing_kmer_length}.{genome_prefix}.{haplotype}.{pairprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bam_merge_pairs.{assembly_stage}.{parameters}.{phasing_kmer_length}.{genome_prefix}.{haplotype}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["two_read_bam_combiner"] ,
        time=parameters["time"]["two_read_bam_combiner"],
        mem=parameters["memory_mb"]["two_read_bam_combiner"] + parameters["memory_mb"]["samtools_sort"]
    threads: parameters["threads"]["two_read_bam_combiner"] + parameters["threads"]["samtools_sort"]
    shell:
        " TMP_PREFIX=`dirname {output.bam}`/{wildcards.pairprefix}; "
        " two_read_bam_combiner.pl {input.forward_bam} {input.reverse_bam} samtools {params.min_mapq} 2>{log.merge} | "
        " samtools view -bS -t {input.reference_fai} - 2>{log.view} | "
        " samtools sort -T ${{TMP_PREFIX}} -m {params.sort_memory}M -@ {params.sort_threads} -o {output.bam} 2>{log.sort}"

rule bam_merge_files:
    input:
        bams=expand(rules.bam_merge_pairs.output.bam, #out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.bwa.filtered.{pairprefix}.bam",
                    allow_missing=True,
                    pairprefix=input_pairprefix_dict["hic"]), #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        reference_fai=rules.ref_faidx.output.fai,
        reference=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta"
    output:
        bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.bwa.bam" # TODO: make temp
    params:
        sort_threads=parameters["threads"]["samtools_sort"]
    log:
        std=output_dict["log"] / "bam_merge_files.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "bam_merge_files.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bam_merge_files.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bam_merge_files.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["samtools_sort"] ,
        time=parameters["time"]["samtools_sort"],
        mem=parameters["memory_mb"]["samtools_sort"]
    threads: parameters["threads"]["samtools_sort"]
    shell:
        " samtools merge -@ {params.sort_threads} -o {output.bam} {input.bams} 1>{log.std} 2>&1"

rule rmdup:
    input:
        bam=rules.bam_merge_files.output.bam
    output:
        bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
        bai=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam.bai",
        dup_stats=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.stats",
        bam_stats=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam.stats",
    log:
        std=output_dict["log"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "rmdup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["rmdup"] ,
        time=parameters["time"]["rmdup"],
        mem=parameters["memory_mb"]["rmdup"]
    threads: parameters["threads"]["rmdup"]
    shell:
        " TMP_DIRECTORY=`dirname {input.bam}`/temp; "
        " picard -Xmx{resources.mem}m MarkDuplicates -I {input} -O {output.bam} " 
        " --REMOVE_DUPLICATES true -M {output.dup_stats}  --TMP_DIR ${{TMP_DIRECTORY}} >{log.std} 2>&1;"
        " samtools index {output.bam}; "
        " get_stats.pl {output.bam} > {output.bam_stats}"

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
        cpus=parameters["threads"]["generate_prescaffolding_agp"] ,
        time=parameters["time"]["generate_prescaffolding_agp"],
        mem=parameters["memory_mb"]["generate_prescaffolding_agp"]
    threads: parameters["threads"]["generate_prescaffolding_agp"]

    shell:
        " ./workflow/scripts/generate_agp_from_fai.py -i {input.reference_fai} -o {output.pre_agp} > {log.std} 2>&1;"

rule yahs_juicer_pre_prescaffolding: #
    input:
        reference_fai=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta.fai",
        bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
        bai=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam.bai",
        agp=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.pre.agp"
    output:
        links_bed=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bed",
        liftover_agp=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.liftover.agp",
        assembly=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.assembly",
        assembly_agp=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.assembly.agp",
        log=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.yahs.juicer_pre.log",
    log:
        #std=output_dict["log"]  / "yahs_juicer_pre.{prev_stage_parameters}..yahs_{hic_scaffolding_parameters}.{genome_prefix}.{haplotype}.log",
        mv=output_dict["log"]  / "yahs_juicer_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.mv.log",
        cluster_log=output_dict["cluster_log"] / "yahs_juicer_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "yahs_juicer_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "yahs_juicer_pre_prescaffolding.{assembly_stage}.{parameters}.{haplotype}.{phasing_kmer_length}.{genome_prefix}.benchmark.txt"
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
        " juicer pre -a -o ${{OUTPUT_PREFIX}} {input.bam} {input.agp} {input.reference_fai} > {output.log} 2>&1;"
        " mv ${{OUTPUT_PREFIX}}.txt {output.links_bed} > {log.mv} 2>&1; "

rule juicer_tools_pre_prescaffolding: #
    input:
        yahs_juicer_pre_log=rules.yahs_juicer_pre_prescaffolding.output.log,
        yahs_juicer_pre_bed=rules.yahs_juicer_pre_prescaffolding.output.links_bed
    output:
        hic=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.hic",
    log:
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
        cpus=parameters["threads"]["juicer_tools_pre"] ,
        time=parameters["time"]["juicer_tools_pre"],
        mem=parameters["memory_mb"]["juicer_tools_pre"]
    threads: parameters["threads"]["juicer_tools_pre"]

    shell: # juicer_tools elder than 1.9.9 seems to be incompartible with yahs
        " java -jar -Xmx{resources.mem}m workflow/external_tools/juicer/juicer_tools.1.9.9_jcuda.0.8.jar pre " #--threads {threads}  
        " {input.yahs_juicer_pre_bed} {output.hic} <(cat {input.yahs_juicer_pre_log} 2>{log.cat} | "
        " grep PRE_C_SIZE 2>{log.grep} | awk '{{print $2\" \"$3}}' 2>{log.awk}) > {log.juicer} 2>&1; "
