

rule bwa_map: #
    input:
        index=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta.ann",
        reference=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",

        forward_fastq=lambda wildcards: output_dict["data"] / "fastq/hic/raw/{0}{1}{2}".format(wildcards.pairprefix,
                                                                                               input_forward_suffix_dict["hic"] if wildcards.phasing_kmer_length == "NA" else "_1",
                                                                                               config["fastq_extension"]) if wildcards.phasing_kmer_length == "NA" else \
                                out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}{6}".format(config["phasing_stage"], #wildcards.assembly_stage,
                                                                                            detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
                                                                                            wildcards.haplotype,
                                                                                            wildcards.phasing_kmer_length,
                                                                                            wildcards.pairprefix,
                                                                                            input_forward_suffix_dict["hic"] if wildcards.phasing_kmer_length == "NA" else "_1",
                                                                                            config["fastq_extension"]),
        reverse_fastq=lambda wildcards: output_dict["data"] / "fastq/hic/raw/{0}{1}{2}".format(wildcards.pairprefix,
                                                                                               input_reverse_suffix_dict["hic"] if wildcards.phasing_kmer_length == "NA" else "_2",
                                                                                               config["fastq_extension"]) if wildcards.phasing_kmer_length == "NA" else \
                                out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}{6}".format(config["phasing_stage"], #wildcards.assembly_stage,
                                                                                            detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
                                                                                            wildcards.haplotype,
                                                                                            wildcards.phasing_kmer_length,
                                                                                            wildcards.pairprefix,
                                                                                            input_reverse_suffix_dict["hic"] if wildcards.phasing_kmer_length == "NA" else "_2",
                                                                                            config["fastq_extension"]),
    output:
        #bam=out_dir_path  / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{fileprefix}.bwa.bam"
        bam=temp(out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{pairprefix}.bwa.bam")
    params:
        id="{0}_hic".format(config["genome_prefix"]),
        bwa_tool=config["bwa_tool"]
    log:
        map=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.map.log",
        sort=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.sort.log",
        #filter=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.filter.log",
        cluster_log=output_dict["cluster_log"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("bwa_map"),
        cpus=parameters["threads"]["bwa_map"] ,
        time=parameters["time"]["bwa_map"],
        mem=parameters["memory_mb"]["bwa_map"]
    threads: parameters["threads"]["bwa_map"]
    shell:
        " {params.bwa_tool} mem -SP5M -t {threads} -R  \'@RG\\tID:{params.id}\\tPU:x\\tSM:{params.id}\\tPL:illumina\\tLB:x\' "
        " {input.reference} {input.forward_fastq} {input.reverse_fastq} 2>{log.map} | samtools view -Sb - > {output.bam} 2>{log.sort} "

rule pairtools_parse:
    input:
        bam=rules.bwa_map.output.bam,
        len_file=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len"
    output:
        pairsam_gz=temp(out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{pairprefix}.bwa.pairsam.gz")
    params:
        min_mapping_quality=lambda wildcards: parse_option("min_mapping_quality", parameters["tool_options"]["pairtools_parse"], " --min-mapq "),
        max_interalign_gap=lambda wildcards: parse_option("max_interalign_gap", parameters["tool_options"]["pairtools_parse"], " --max-inter-align-gap ")
    log:
        std=output_dict["log"] / "pairtools_parse.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.log",
        cluster_log=output_dict["cluster_log"] / "pairtools_parse.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pairtools_parse.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pairtools_parse.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("pairtools_parse"),
        cpus=parameters["threads"]["pairtools_parse"] ,
        time=parameters["time"]["pairtools_parse"],
        mem=parameters["memory_mb"]["pairtools_parse"]
    threads: parameters["threads"]["pairtools_parse"]
    shell:
        " samtools view -h {input.bam} | "
        " pairtools parse {params.min_mapping_quality} --walks-policy 5unique {params.max_interalign_gap} "
        " --nproc-in {threads} --nproc-out {threads} --chroms-path {input.len_file} -o {output.pairsam_gz} > {log.std} 2>&1; "

rule pairtools_sort:
    input:
        pairsam_gz=rules.pairtools_parse.output.pairsam_gz
    output:
        sorted_pairsam_gz=temp(out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{pairprefix}.bwa.sorted.pairsam.gz")
    log:
        std=output_dict["log"] / "pairtools_sort.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.log",
        cluster_log=output_dict["cluster_log"] / "pairtools_sort.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pairtools_sort.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pairtools_sort.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("pairtools_sort"),
        cpus=parameters["threads"]["pairtools_sort"] ,
        time=parameters["time"]["pairtools_sort"],
        mem=parameters["memory_mb"]["pairtools_sort"]
    threads: parameters["threads"]["pairtools_sort"]
    shell:
        " TMP_DIR=`dirname {output.sorted_pairsam_gz}`/{wildcards.pairprefix}_tmp; "
        " mkdir -p ${{TMP_DIR}}; "
        " pairtools sort --nproc {threads} --memory {resources.mem}M --tmpdir=${{TMP_DIR}} "
        " -o {output.sorted_pairsam_gz} {input.pairsam_gz}  > {log.std} 2>&1; "

rule pairtools_merge:
    input:
        pairsam_gzs=expand(rules.pairtools_parse.output.pairsam_gz,
                           pairprefix=input_pairprefix_dict["hic"],
                           allow_missing=True)
    output:
        merged_pairsam_gz=temp(out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.bwa.merged.pairsam.gz")
    params:
        memory=int(0.5 * parameters["memory_mb"]["pairtools_merge"])
    log:
        std=output_dict["log"] / "pairtools_merge.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "pairtools_merge.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pairtools_merge.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pairtools_merge.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("pairtools_merge"),
        cpus=parameters["threads"]["pairtools_merge"] ,
        time=parameters["time"]["pairtools_merge"],
        mem=parameters["memory_mb"]["pairtools_merge"]
    threads: parameters["threads"]["pairtools_merge"]
    shell:
        " TMP_DIR=`dirname {output.merged_pairsam_gz}`/merged_tmp; "
        " mkdir -p ${{TMP_DIR}}; "
        " pairtools merge --nproc {threads} --max-nmerge 8 --memory {params.memory}M --tmpdir=${{TMP_DIR}} "
        " -o {output.merged_pairsam_gz} {input.pairsam_gzs}  > {log.std} 2>&1; "

rule pairtools_dedup:
    input:
        merged_pairsam_gz=rules.pairtools_merge.output.merged_pairsam_gz
    output:
        dedup_pairsam_gz=temp(out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.bwa.rmdup.pairsam.gz"),
        dedup_pairsam_stats=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.bwa.rmdup.pairsam.stats",
        dedup_pairsam_summary=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.bwa.rmdup.pairsam.summary"
    log:
        std=output_dict["log"] / "pairtools_dedup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.log",
        summary=output_dict["log"] / "pairtools_dedup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.summary.log",
        cluster_log=output_dict["cluster_log"] / "pairtools_dedup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pairtools_dedup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pairtools_dedup.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
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
        sorted_dedup_bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
        pairs=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.pairs.gz"
    params:
        sort_threads=parameters["threads"]["samtools_sort"],
        sort_per_thread=parameters["memory_mb"]["samtools_sort"],
        split_threads=parameters["threads"]["pairtools_split"],
        sort_mem=parameters["memory_mb"]["samtools_sort"] * parameters["threads"]["samtools_sort"],
    log:
        split=output_dict["log"] / "pairtools_split.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.split.log",
        sort=output_dict["log"] / "pairtools_split.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.sort.log",
        sort_pairs=output_dict["log"] / "pairtools_split.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.sort_pairs.log",
        cluster_log=output_dict["cluster_log"] / "pairtools_split.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pairtools_split.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pairtools_split.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
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
        " --output-pairs ${{UNSORTED_PAIRS}} --output-sam - {input.dedup_pairsam_gz} 2>{log.split} | "
        " samtools sort -@ {params.sort_threads} -m {params.sort_per_thread}M -T ${{TMP_PREFIX}} -o {output.sorted_dedup_bam} > {log.sort} 2>&1; "
        " pairtools sort --nproc-in {threads} --nproc-out {threads} --memory {params.sort_mem}M --tmpdir ${{TMP_DIR}} "
        " -o {output.pairs} ${{UNSORTED_PAIRS}} > {log.sort_pairs} 2>&1 ; "

rule pairtools_index_pairs:
    input:
        pairs=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.pairs.gz"
    output:
        index=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.pairs.gz.px2"
        #sorted_dedup_bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam",
    log:
        std=output_dict["log"] / "pairtools_index_pairs.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.std.log",
        cluster_log=output_dict["cluster_log"] / "pairtools_index_pairs.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pairtools_index_pairs.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pairtools_index_pairs.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("pairtools_index_pairs"),
        cpus=parameters["threads"]["pairtools_index"],
        time=parameters["time"]["pairtools_index"],
        mem=parameters["memory_mb"]["pairtools_index"],
    threads: parameters["threads"]["pairtools_index"]
    shell:
        " pairix -f -p pairs {input.pairs} > {log.std} 2>&1; "