

rule nextpolish2: #
    input:
        draft_assembly=lambda wildcards: config["out_dir"] / "{2}/{0}/{1}.{2}.{3}.fasta".format(get_prev_stage_parameters(wildcards.parameters),
                                                                                                 wildcards.genome_prefix,
                                                                                                 stage_dict["polishing"].prev_stage,
                                                                                                 wildcards.haplotype),
        winnowmap_bam=lambda wildcards: config["out_dir"] / "{2}/{0}/alignment/hifi/winnowmap/{1}.{2}.{3}.k15.hifi.winnowmap.bam".format(get_prev_stage_parameters(wildcards.parameters),
                                                                                                 wildcards.genome_prefix,
                                                                                                 stage_dict["polishing"].prev_stage,
                                                                                                 wildcards.haplotype),
        winnowmap_bam_index=lambda wildcards: config["out_dir"] / "{2}/{0}/alignment/hifi/winnowmap/{1}.{2}.{3}.k15.hifi.winnowmap.bam.csi".format(get_prev_stage_parameters(wildcards.parameters),
                                                                                                 wildcards.genome_prefix,
                                                                                                 stage_dict["polishing"].prev_stage,
                                                                                                 wildcards.haplotype),
        yak_db_k31=config["out_dir"] / "kmer/illumina/final/illumina.final.31.yak",
        yak_db_k21=config["out_dir"] / "kmer/illumina/final/illumina.final.21.yak",
        log_dir=ancient(config["out_dir"] / "polishing/{parameters}/log/"),
    output:
        polished_assembly=config["out_dir"] / "polishing/{parameters, [^/]*nextpolish2[^/]*}/{genome_prefix}.polishing.{haplotype, hap[^/]+}.fasta"
    params:
        sort_threads=parameters["threads"]["samtools_sort"],
        minimap_threads=parameters["threads"]["minimap2"],
        per_thread_sort_mem=parameters["memory_mb"]["samtools_sort"],
    log:
        std=config["out_dir"] / "polishing/{parameters}/log/nextpolish2.hifi.{parameters}.{genome_prefix}.polishing.{haplotype}.map.log",
        cluster_log=config["out_dir"] / "polishing/{parameters}/log/nextpolish2.hifi.{parameters}.{genome_prefix}.polishing.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "polishing/{parameters}/log/nextpolish2.hifi.{parameters}.{genome_prefix}.polishing.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "polishing/{parameters}/log/nextpolish2.hifi.{parameters}.{genome_prefix}.polishing.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["nextpolish2"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["nextpolish2"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("nextpolish2"),
        cpus=parameters["threads"]["nextpolish2"],
        time=parameters["time"]["nextpolish2"],
        mem=parameters["memory_mb"]["nextpolish2"]
    threads: parameters["threads"]["nextpolish2"]
    shell:
        " nextPolish2  -t {threads} {input.winnowmap_bam} {input.draft_assembly} "
        "     {input.yak_db_k21} {input.yak_db_k31} > {output.polished_assembly} 2>{log.std};  "
