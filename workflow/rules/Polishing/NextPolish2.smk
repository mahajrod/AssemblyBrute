

rule nextpolish2: #
    input:
        draft_assembly=lambda wildcards: out_dir_path / "{2}/{0}/{1}.{2}.{3}.fasta".format(wildcards.prev_stage_parameters,
                                                                                                 wildcards.genome_prefix,
                                                                                                 stage_dict["polishing"]["parameters"][wildcards.prev_stage_parameters + "..nextpolish2_" + wildcards.polishing_parameters]["prev_stage"],
                                                                                                 wildcards.haplotype) ,
        winnowmap_bam=lambda wildcards: out_dir_path / "{2}/{0}/alignment/hifi/winnowmap/{1}.{2}.{3}.k15.hifi.winnowmap.bam".format(wildcards.prev_stage_parameters,
                                                                                                 wildcards.genome_prefix,
                                                                                                 stage_dict["polishing"]["parameters"][wildcards.prev_stage_parameters + "..nextpolish2_" + wildcards.polishing_parameters]["prev_stage"],
                                                                                                 wildcards.haplotype) ,
        winnowmap_bam_index=lambda wildcards: out_dir_path / "{2}/{0}/alignment/hifi/winnowmap/{1}.{2}.{3}.k15.hifi.winnowmap.bam.csi".format(wildcards.prev_stage_parameters,
                                                                                                 wildcards.genome_prefix,
                                                                                                 stage_dict["polishing"]["parameters"][wildcards.prev_stage_parameters + "..nextpolish2_" + wildcards.polishing_parameters]["prev_stage"],
                                                                                                 wildcards.haplotype) ,
        yak_db_k31=output_dict["kmer"] / "illumina/{0}/illumina.{0}.31.yak".format("filtered" if "illumina" in config["filtered_data"] else "raw"),
        yak_db_k21=output_dict["kmer"] / "illumina/{0}/illumina.{0}.21.yak".format("filtered" if "illumina" in config["filtered_data"] else "raw"),
        log_dir=out_dir_path / "polishing/{prev_stage_parameters}..nextpolish2_{polishing_parameters}/log/",
        benchmark_dir=out_dir_path / "polishing/{prev_stage_parameters}..nextpolish2_{polishing_parameters}/benchmark/"
    output:
        #bam=out_dir_path  / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{fileprefix}.bwa.bam"
        polished_assembly=out_dir_path / "polishing/{prev_stage_parameters, [^/]+}..nextpolish2_{polishing_parameters, [^/]+}/{genome_prefix, [^/]+}.polishing.{haplotype, hap[^/]+}.fasta"
    params:
        sort_threads=parameters["threads"]["samtools_sort"],
        minimap_threads=parameters["threads"]["minimap2"],
        per_thread_sort_mem=parameters["memory_mb"]["samtools_sort"],
    log:
        std=out_dir_path / "polishing/{prev_stage_parameters}..nextpolish2_{polishing_parameters}/log/nextpolish2.hifi.{genome_prefix}.polishing.{haplotype}.map.log",
        cluster_log=out_dir_path / "polishing/{prev_stage_parameters}..nextpolish2_{polishing_parameters}/log/nextpolish2.hifi.{genome_prefix}.polishing.{haplotype}.cluster.log",
        cluster_err=out_dir_path / "polishing/{prev_stage_parameters}..nextpolish2_{polishing_parameters}/log/nextpolish2.hifi.{genome_prefix}.polishing.{haplotype}.cluster.err"
    benchmark:
        out_dir_path / "polishing/{prev_stage_parameters}..nextpolish2_{polishing_parameters}/benchmark/nextpolish2.hifi.{genome_prefix}.polishing.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("nextpolish2"),
        cpus=parameters["threads"]["nextpolish2"],
        time=parameters["time"]["nextpolish2"],
        mem=parameters["memory_mb"]["nextpolish2"]
    threads: parameters["threads"]["nextpolish2"]
    shell:
        " nextPolish2 -r -t {threads} {input.winnowmap_bam} {input.draft_assembly} "
        "             {input.yak_db_k21} {input.yak_db_k31} > asm.np2.fa 2>{log.std};  "


