ruleorder: arima_two_read_bam_combiner > arima_filter_five_end

rule arima_bwa_map: #
    input:
        index=config["out_dir"] / "hic_alignment/{parameters}/{genome_prefix}.hic_alignment.{haplotype}.fasta.ann",
        reference=config["out_dir"] / "hic_alignment/{parameters}/{genome_prefix}.hic_alignment.{haplotype}.fasta",
        fastq=lambda wildcards: config["out_dir"] / "data/hic/final/{0}{1}".format(wildcards.fileprefix,
                                                                                    config["data"]["hic"]["conv_ext"]) if wildcards.phasing_kmer_length == "NA" else \
                                config["out_dir"] / "{0}/{1}/{6}.{0}.{2}/reads/hic/{3}/{4}{5}".format(config["phasing_stage"],
                                                                                                      detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."),
                                                                                                      wildcards.haplotype,
                                                                                                      wildcards.phasing_kmer_length,
                                                                                                      wildcards.fileprefix,
                                                                                                      config["data"]["hic"]["conv_ext"],
                                                                                                      config["genome_prefix"]
                                                                                                      )
    output:
        raw_bam=temp(config["out_dir"]  / "hic_alignment/{parameters, .*arima.*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{fileprefix}.bwa.raw.bam")
    log:
        fastx=config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.fastx.log",
        map=config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.map.log",
        sort=config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.sort.log",
        filter=config["out_dir"]  / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.filter.log",
        cluster_log=config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{fileprefix}.cluster.err"
    benchmark:
        config["out_dir"]  / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("bwa_map"),
        cpus=get_threads(parameters["threads"]["bwa_map_arima"], "cpu"),
        time=parameters["time"]["bwa_map"],
        mem=parameters["memory_mb"]["bwa_map"],
    threads: parameters["threads"]["bwa_map_arima"]
    shell:
        " bwa-mem2 mem -SP5M -t {threads} -R  \'@RG\\tID:{wildcards.genome_prefix}_hic\\tPU:x\\tSM:{wildcards.genome_prefix}_hic\\tPL:illumina\\tLB:x\' "
        "     {input.reference} <(zcat {input.fastq}  2>{log.fastx}) 2>{log.map} | samtools view -Sb - > {output.raw_bam} 2>{log.sort} "

rule arima_filter_five_end: #
    input:
        raw_bam=config["out_dir"] / "hic_alignment/{parameters}/{genome_prefix}.hic_alignment.{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{fileprefix}.bwa.raw.bam",
        stats=config["out_dir"]  / "hic_alignment/{parameters}/{genome_prefix}.hic_alignment.{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{fileprefix}.bwa.raw.bam.general_stats"
    output:
        bam=temp(config["out_dir"]  / "hic_alignment/{parameters, .*arima.*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{fileprefix}.arima.bwa.bam")
    params:
        id="{0}_hic".format(config["genome_prefix"]),
    log:
        view=config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.view.log",
        filter=config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.filter.log",
        view2=config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.view2.log",
        cluster_log=config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{fileprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("bwa_map"),
        cpus=parameters["threads"]["arima_filter_five_end"] ,
        time=parameters["time"]["arima_filter_five_end"],
        mem=parameters["memory_mb"]["arima_filter_five_end"],
    threads: parameters["threads"]["arima_filter_five_end"]
    shell:
        " samtools view -h {input.raw_bam} 2>{log.view} | "
        "     workflow/external_tools/arima_mapping_pipeline/filter_five_end.pl 2>{log.filter} | "
        "     samtools view -Sb - > {output.bam} 2>{log.view2} "


rule arima_two_read_bam_combiner:
    input:
        forward_bam=lambda wildcards: expand(rules.arima_filter_five_end.output.bam,
                                             fileprefix=[wildcards.pairprefix + config["data"]["hic"]["conv_fwd_sfx"]],
                                             allow_missing=True),
        reverse_bam=lambda wildcards: expand(rules.arima_filter_five_end.output.bam,
                                             fileprefix=[wildcards.pairprefix + config["data"]["hic"]["conv_rev_sfx"]],
                                            allow_missing=True),
        reference_fai=config["out_dir"] / "hic_alignment/{parameters}/{genome_prefix}.hic_alignment.{haplotype}.fasta.fai"
    output:
        bam=temp(config["out_dir"] / "hic_alignment/{parameters, .*arima.*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{pairprefix}.bwa.bam"),
    params:
        min_mapq=parameters["tool_options"]["two_read_bam_combiner"]["mapq"],
        sort_threads=parameters["threads"]["samtools_sort"],
        sort_memory=parameters["memory_mb"]["samtools_sort"],
    log:
        merge=config["out_dir"] / "log/bam_merge_pairs.{parameters}.{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{pairprefix}.merge.log",
        view=config["out_dir"] / "log/bam_merge_pairs.{parameters}.{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{pairprefix}.view.log",
        sort=config["out_dir"] / "log/bam_merge_pairs.{parameters}.{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{pairprefix}.sort.log",
        cluster_log=config["out_dir"] / "log/bam_merge_pairs.{parameters}.{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{pairprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/bam_merge_pairs.{parameters}.{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{pairprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/bam_merge_pairs.{parameters}.{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["arima_two_read_bam_combiner"] ,
        node_options=parse_node_list("bwa_merge_pairs"),
        time=parameters["time"]["arima_two_read_bam_combiner"],
        mem=parameters["memory_mb"]["arima_two_read_bam_combiner"] + parameters["memory_mb"]["samtools_sort"] * parameters["threads"]["samtools_sort"] + 30000,
        queue=config["queue"]["cpu"]["name"]
    threads: parameters["threads"]["arima_two_read_bam_combiner"] + parameters["threads"]["samtools_sort"]
    shell:
        " TMP_PREFIX=`dirname {output.bam}`/{wildcards.pairprefix}; "
        " workflow/external_tools/arima_mapping_pipeline/two_read_bam_combiner.pl {input.forward_bam} "
        "     {input.reverse_bam} samtools {params.min_mapq} 2>{log.merge} | "
        "     samtools view -bS -t {input.reference_fai} - 2>{log.view} | "
        "     samtools sort -T ${{TMP_PREFIX}} -m {params.sort_memory}M -@ {params.sort_threads} -o {output.bam} 2>{log.sort}"
