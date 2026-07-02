ruleorder: pairtools_split > rmdup

rule pairtools_parse:
    input:
        bam=rules.bwa_map.output.bam,
        len_file=config["out_dir"]  / "hic_alignment/{parameters}/{genome_prefix}.hic_alignment.{haplotype}.len"
    output:
        pairsam_gz=temp(config["out_dir"] / "hic_alignment/{parameters, .*pairtools.*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{pairprefix}.bwa.pairsam.gz")
    params:
        min_mapping_quality=lambda wildcards: parse_option("min_mapping_quality", parameters["tool_options"]["pairtools_parse"], " --min-mapq "),
        max_interalign_gap=lambda wildcards: parse_option("max_interalign_gap", parameters["tool_options"]["pairtools_parse"], " --max-inter-align-gap ")
    log:
        std=config["out_dir"] / "log/pairtools_parse.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.log",
        cluster_log=config["out_dir"] / "log/pairtools_parse.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/pairtools_parse.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/pairtools_parse.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["pairtools"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["pairtools"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("pairtools_parse"),
        cpus=parameters["threads"]["pairtools_parse"] ,
        time=parameters["time"]["pairtools_parse"],
        mem=parameters["memory_mb"]["pairtools_parse"]
    threads: parameters["threads"]["pairtools_parse"]
    shell:
        " samtools view -h {input.bam} | "
        " pairtools parse {params.min_mapping_quality} --walks-policy 5unique {params.max_interalign_gap} "
        "     --nproc-in {threads} --nproc-out {threads} --chroms-path {input.len_file} -o {output.pairsam_gz} > {log.std} 2>&1; "

rule pairtools_sort:
    input:
        pairsam_gz=rules.pairtools_parse.output.pairsam_gz
    output:
        sorted_pairsam_gz=temp(config["out_dir"] / "hic_alignment/{parameters, .*pairtools.*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{pairprefix}.bwa.sorted.pairsam.gz")
    log:
        std=config["out_dir"] / "log/pairtools_sort.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.log",
        cluster_log=config["out_dir"] / "log/pairtools_sort.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/pairtools_sort.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/pairtools_sort.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["pairtools"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["pairtools"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("pairtools_sort"),
        cpus=parameters["threads"]["pairtools_sort"] ,
        time=parameters["time"]["pairtools_sort"],
        mem=parameters["memory_mb"]["pairtools_sort"]
    threads: parameters["threads"]["pairtools_sort"]
    shell:
        " TMP_DIR=`dirname {output.sorted_pairsam_gz}`/{wildcards.pairprefix}_tmp; "
        " mkdir -p ${{TMP_DIR}}; "
        " pairtools sort --nproc {threads} --memory {resources.mem}M --tmpdir=${{TMP_DIR}} "
        "     -o {output.sorted_pairsam_gz} {input.pairsam_gz}  > {log.std} 2>&1; "

rule pairtools_merge:
    input:
        pairsam_gzs=expand(rules.pairtools_parse.output.pairsam_gz,
                           pairprefix=config["data"]["hic"]["pair_prefix_list"],
                           allow_missing=True)
    output:
        merged_pairsam_gz=temp(config["out_dir"] / "hic_alignment/{parameters, .*pairtools.*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.bwa.merged.pairsam.gz")
    params:
        memory=int(0.5 * parameters["memory_mb"]["pairtools_merge"])
    log:
        std=config["out_dir"] / "log/pairtools_merge.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.log",
        cluster_log=config["out_dir"] / "log/pairtools_merge.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/pairtools_merge.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/pairtools_merge.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["pairtools"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["pairtools"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("pairtools_merge"),
        cpus=parameters["threads"]["pairtools_merge"] ,
        time=parameters["time"]["pairtools_merge"],
        mem=parameters["memory_mb"]["pairtools_merge"]
    threads: parameters["threads"]["pairtools_merge"]
    shell:
        " TMP_DIR=`dirname {output.merged_pairsam_gz}`/merged_tmp; "
        " mkdir -p ${{TMP_DIR}}; "
        " pairtools merge --nproc {threads} --max-nmerge 8 --memory {params.memory}M --tmpdir=${{TMP_DIR}} "
        "     -o {output.merged_pairsam_gz} {input.pairsam_gzs}  > {log.std} 2>&1; "

rule pairtools_dedup:
    input:
        merged_pairsam_gz=rules.pairtools_merge.output.merged_pairsam_gz
    output:
        dedup_pairsam_gz=temp(config["out_dir"] / "hic_alignment/{parameters, .*pairtools.*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.bwa.rmdup.pairsam.gz"),
        dedup_pairsam_stats=config["out_dir"] / "hic_alignment/{parameters, .*pairtools.*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.bwa.rmdup.pairsam.stats",
        dedup_pairsam_summary=config["out_dir"] / "hic_alignment/{parameters, .*pairtools.*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}{phasing_kmer_length}..bwa.rmdup.pairsam.summary"
    log:
        std=config["out_dir"] / "log/pairtools_dedup.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.log",
        summary=config["out_dir"] / "log/pairtools_dedup.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.summary.log",
        cluster_log=config["out_dir"] / "log/pairtools_dedup.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/pairtools_dedup.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/pairtools_dedup.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["pairtools"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["pairtools"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("pairtools_dedup"),
        cpus=parameters["threads"]["pairtools_dedup"] ,
        time=parameters["time"]["pairtools_dedup"],
        mem=parameters["memory_mb"]["pairtools_dedup"]
    threads: parameters["threads"]["pairtools_dedup"]
    shell:
        " pairtools dedup --nproc-in {threads} --nproc-out {threads} --mark-dups --output-stats {output.dedup_pairsam_stats} "
        " --output {output.dedup_pairsam_gz} {input.merged_pairsam_gz} > {log.std} 2>&1; "
        " workflow/external_tools/Omni-C/get_qc.py -p {output.dedup_pairsam_stats} > {output.dedup_pairsam_summary} 2>{log.summary}; "

rule pairtools_split:
    input:
        dedup_pairsam_gz=rules.pairtools_dedup.output.dedup_pairsam_gz
    output:
        sorted_dedup_bam=config["out_dir"] / "hic_alignment/{parameters, .*pairtools.*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.rmdup.bam",
        pairs=config["out_dir"] / "hic_alignment/{parameters, .*pairtools.*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.rmdup.pairs.gz"
    params:
        sort_threads=parameters["threads"]["samtools_sort"],
        sort_per_thread=parameters["memory_mb"]["samtools_sort"],
        split_threads=parameters["threads"]["pairtools_split"],
        sort_mem=parameters["memory_mb"]["samtools_sort"] * parameters["threads"]["samtools_sort"],
    log:
        split=config["out_dir"] / "log/pairtools_split.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.split.log",
        sort=config["out_dir"] / "log/pairtools_split.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.sort.log",
        sort_pairs=config["out_dir"] / "log/pairtools_split.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.sort_pairs.log",
        cluster_log=config["out_dir"] / "log/pairtools_split.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/pairtools_split.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/pairtools_split.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["pairtools"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["pairtools"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("pairtools_split"),
        cpus=parameters["threads"]["pairtools_split"] + parameters["threads"]["samtools_sort"],
        time=parameters["time"]["pairtools_split"],
        mem=parameters["memory_mb"]["pairtools_split"] + parameters["memory_mb"]["samtools_sort"] * parameters["threads"]["samtools_sort"]
    threads: parameters["threads"]["pairtools_split"]
    shell:
        " TMP_DIR=`dirname {output.pairs}`; "
        " TMP_PREFIX=${{TMP_DIR}}/pairtools_samtools_sort_tmp; "
        " UNSORTED_PAIRS={output.pairs}; "
        " UNSORTED_PAIRS=${{UNSORTED_PAIRS%.gz}}; "
        " pairtools split --nproc-in {params.split_threads} --nproc-out {params.split_threads} "
        "     --output-pairs ${{UNSORTED_PAIRS}} --output-sam - {input.dedup_pairsam_gz} 2>{log.split} | "
        "     samtools sort -@ {params.sort_threads} -m {params.sort_per_thread}M -T ${{TMP_PREFIX}} -o {output.sorted_dedup_bam} > {log.sort} 2>&1; "
        " pairtools sort --nproc-in {threads} --nproc-out {threads} --memory {params.sort_mem}M --tmpdir ${{TMP_DIR}} "
        "     -o {output.pairs} ${{UNSORTED_PAIRS}} > {log.sort_pairs} 2>&1 ; "

rule pairtools_index_pairs:
    input:
        pairs=config["out_dir"] / "hic_alignment/{parameters}/{genome_prefix}.hic_alignment.{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.rmdup.pairs.gz"
    output:
        index=config["out_dir"] / "hic_alignment/{parameters, .*pairtools.*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.rmdup.pairs.gz.px2"
    log:
        std=config["out_dir"] / "log/pairtools_index_pairs.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.std.log",
        cluster_log=config["out_dir"] / "log/pairtools_index_pairs.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/pairtools_index_pairs.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/pairtools_index_pairs.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["pairtools"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["pairtools"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("pairtools_index_pairs"),
        cpus=parameters["threads"]["pairtools_index"],
        time=parameters["time"]["pairtools_index"],
        mem=parameters["memory_mb"]["pairtools_index"],
    threads: parameters["threads"]["pairtools_index"]
    shell:
        " pairix -f -p pairs {input.pairs} > {log.std} 2>&1; "
