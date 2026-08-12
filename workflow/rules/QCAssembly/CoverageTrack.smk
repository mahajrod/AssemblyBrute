ruleorder: get_track_stats > create_coverage_table
ruleorder: create_bedgraph_from_coverage_table > create_bedgraph_track
localrules: gather_coverage_track_files
if "purge_dups" in config["stage_list"]:
    ruleorder: minimap2_cov > minimap2_purge_dups_reads


rule minimap2_cov:
    input:
        fastq=lambda wildcards: expand(config["out_dir"] / ("data/%s/final/{fileprefix}%s" % (wildcards.longread_datatype,
                                                                                              config["data"][wildcards.longread_datatype]["conv_ext"])),
                     fileprefix=config["data"][wildcards.longread_datatype]["conv_file_prefix_list"],
                     allow_missing=True),
        reference="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        bam="{fasta_dir}/{fasta_prefix}/alignment/{fasta_prefix}.{longread_datatype}.bam"

    params:
        index_size=lambda wildcards: parse_option("index_size", parameters["tool_options"]["minimap2"][wildcards.longread_datatype], " -I "),
        alignment_scheme=lambda wildcards: parse_option("alignment_scheme", parameters["tool_options"]["minimap2"][wildcards.longread_datatype], " -x "),
        sort_threads=parameters["threads"]["samtools_sort"],
        minimap_threads=parameters["threads"]["minimap2"],
        per_thread_sort_mem=parameters["memory_mb"]["samtools_sort"],
    log:
        minimap2="{fasta_dir}/log/minimap2_cov.{fasta_prefix}.{longread_datatype}.minimap2.log",
        sort="{fasta_dir}/log/minimap2_cov.{fasta_prefix}.{longread_datatype}.sort.log",
        index="{fasta_dir}/log/minimap2_cov.{fasta_prefix}.{longread_datatype}.index.log",
        cluster_log="{fasta_dir}/log/minimap2_cov.{fasta_prefix}.{longread_datatype}.cluster.log",
        cluster_err="{fasta_dir}/log/minimap2_cov.{fasta_prefix}.{longread_datatype}.cluster.err"
    benchmark:
        "{fasta_dir}/log/minimap2_cov.{fasta_prefix}.{longread_datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("minimap2_cov"),
        cpus=parameters["threads"]["minimap2"] + parameters["threads"]["samtools_sort"],
        time=parameters["time"]["minimap2"],
        mem=parameters["memory_mb"]["minimap2"] + (parameters["memory_mb"]["samtools_sort"] * parameters["threads"]["samtools_sort"])
    threads: parameters["threads"]["minimap2"] + parameters["threads"]["samtools_sort"]

    shell:
        " TMPDIR=`dirname {output.bam}`; "
        " minimap2 -L {params.alignment_scheme} {params.index_size} -a -t {params.minimap_threads}  {input.reference}  "
        "     {input.fastq} 2>{log.minimap2} |  "
        "     samtools sort -T ${{TMPDIR}} -@ {params.sort_threads} "
        "     -m {params.per_thread_sort_mem}M -o {output.bam} 2>{log.sort};"

use rule minimap2_cov as minimap2_cov_track_data with:
    input:
        fastq=lambda wildcards: expand(config["out_dir"] / ("track_data/%s/%s/final/{fileprefix}%s" % (wildcards.longread_datatype, wildcards.track_name,
                                                                                                       config["track_data"][wildcards.longread_datatype][wildcards.track_name]["conv_ext"])),
                     fileprefix=config["track_data"][wildcards.longread_datatype][wildcards.track_name]["conv_file_prefix_list"],
                     allow_missing=True),
        reference="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        bam="{fasta_dir}/{fasta_prefix}/alignment/{fasta_prefix}.{longread_datatype}_{track_name}.bam"
    log:
        minimap2="{fasta_dir}/log/minimap2_cov_track_data.{fasta_prefix}.{longread_datatype}.{track_name}.minimap2.log",
        sort="{fasta_dir}/log/minimap2_cov_track_data.{fasta_prefix}.{longread_datatype}.{track_name}.sort.log",
        index="{fasta_dir}/log/minimap2_cov_track_data.{fasta_prefix}.{longread_datatype}.{track_name}.index.log",
        cluster_log="{fasta_dir}/log/minimap2_cov_track_data.{fasta_prefix}.{longread_datatype}.{track_name}.cluster.log",
        cluster_err="{fasta_dir}/log/minimap2_cov_track_data.{fasta_prefix}.{longread_datatype}.{track_name}.cluster.err"
    benchmark:
        "{fasta_dir}/log/minimap2_cov_track_data.{fasta_prefix}.{longread_datatype}.{track_name}.benchmark.txt"

rule bwa_cov:
    input:
        forward_fastqs=lambda wildcards: expand(config["out_dir"] / ("data/%s/final/{pairprefix}%s%s" % (wildcards.pe_datatype,
                                                                                                                   config["data"][wildcards.pe_datatype]["conv_fwd_sfx"],
                                                                                                                   config["data"][wildcards.pe_datatype]["conv_ext"])),
                     pairprefix=config["data"][wildcards.pe_datatype]["pair_prefix_list"],
                     allow_missing=True),
        reverse_fastqs=lambda wildcards: expand(config["out_dir"] / ("data/%s/final/{pairprefix}%s%s" % (wildcards.pe_datatype,
                                                                                                                   config["data"][wildcards.pe_datatype]["conv_rev_sfx"],
                                                                                                                   config["data"][wildcards.pe_datatype]["conv_ext"])),
                     pairprefix=config["data"][wildcards.pe_datatype]["pair_prefix_list"],
                     allow_missing=True),
        reference="{fasta_dir}/{fasta_prefix}.fasta",
        reference_index="{fasta_dir}/{fasta_prefix}.fasta.bwt.2bit.64",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        bam="{fasta_dir}/{fasta_prefix}/alignment/{fasta_prefix}.{pe_datatype, illumina}.bam"

    params:
        bwa_threads=parameters["threads"]["bwa_map"],
        sort_threads=parameters["threads"]["samtools_sort"],
        fixmate_threads=parameters["threads"]["samtools_fixmate"],
        markdup_threads=parameters["threads"]["samtools_markdup"],
        per_thread_sort_mem=parameters["memory_mb"]["samtools_sort_per_thread"],
        genome_prefix=config["genome_prefix"]
    log:
        bwa="{fasta_dir}/log/bwa_cov.{fasta_prefix}.{pe_datatype}.bwa.log",
        fixmate="{fasta_dir}/log/bwa_cov.{fasta_prefix}.{pe_datatype}.fixmate.log",
        sort="{fasta_dir}/log/bwa_cov.{fasta_prefix}.{pe_datatype}.sort.log",
        markdup="{fasta_dir}/log/bwa_cov.{fasta_prefix}.{pe_datatype}.markdup.log",
        cluster_log="{fasta_dir}/log/bwa_cov.{fasta_prefix}.{pe_datatype}.cluster.log",
        cluster_err="{fasta_dir}/log/bwa_cov.{fasta_prefix}.{pe_datatype}.cluster.err"
    benchmark:
        "{fasta_dir}/log/bwa_cov.{fasta_prefix}.{pe_datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("bwa_cov"),
        cpus=get_threads(parameters["threads"]["bwa_map"] + parameters["threads"]["samtools_sort"] + parameters["threads"]["samtools_fixmate"] + parameters["threads"]["samtools_markdup"], "cpu"),
        time=parameters["time"]["bwa_map"],
        mem=parameters["memory_mb"]["bwa_map"] + parameters["memory_mb"]["samtools_sort_per_thread"]*parameters["threads"]["samtools_sort"] + parameters["memory_mb"]["samtools_fixmate"] + parameters["memory_mb"]["samtools_markdup"],
    threads: parameters["threads"]["bwa_map"] + parameters["threads"]["samtools_sort"] + parameters["threads"]["samtools_fixmate"] + parameters["threads"]["samtools_markdup"]

    shell:
        " TMP_PREFIX=`dirname {output.bam}`/tmpbam; "
        " bwa-mem2 mem  -t {params.bwa_threads} {input.reference} <(gunzip -c {input.forward_fastqs}) <(gunzip -c {input.reverse_fastqs}) "
        "     -R  \'@RG\\tID:{params.genome_prefix}\\tPU:x\\tSM:{params.genome_prefix}\\tPL:Illumina\\tLB:x\' 2>{log.bwa} | "
        "     samtools fixmate -@ {params.fixmate_threads} -m - -  2>{log.fixmate} | "
        "     samtools sort -T {{TMP_PREFIX}} -@ {params.sort_threads} -m {params.per_thread_sort_mem}M 2>{log.sort} | "
        "     samtools markdup -@ {params.markdup_threads} - {output.bam} 2>{log.markdup}"

use rule bwa_cov as bwa_cov_track_data with:
    input:
        forward_fastqs=lambda wildcards: expand(config["out_dir"] / ("track_data/%s/%s/final/{pairprefix}%s%s" % (wildcards.pe_datatype, wildcards.track_name,
                                                                                                                   config["track_data"][wildcards.pe_datatype][wildcards.track_name]["conv_fwd_sfx"],
                                                                                                                   config["track_data"][wildcards.pe_datatype][wildcards.track_name]["conv_ext"])),
                     pairprefix=config["data"][wildcards.pe_datatype]["pair_prefix_list"],
                     allow_missing=True),
        reverse_fastqs=lambda wildcards: expand(config["out_dir"] / ("track_data/%s/%s/final/{pairprefix}%s%s" % (wildcards.pe_datatype, wildcards.track_name,
                                                                                                               config["track_data"][wildcards.pe_datatype][wildcards.track_name]["conv_rev_sfx"],
                                                                                                               config["track_data"][wildcards.pe_datatype][wildcards.track_name]["conv_ext"])),
                     pairprefix=config["track_data"][wildcards.pe_datatype][wildcards.track_name]["pair_prefix_list"],
                     allow_missing=True),
        reference="{fasta_dir}/{fasta_prefix}.fasta",
        reference_index="{fasta_dir}/{fasta_prefix}.fasta.bwt.2bit.64",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        bam="{fasta_dir}/{fasta_prefix}/alignment/{fasta_prefix}.{pe_datatype, illumina}_{track_name}.bam"
    log:
        bwa="{fasta_dir}/log/bwa_cov_track_data.{fasta_prefix}.{pe_datatype}.{track_name}.bwa.log",
        fixmate="{fasta_dir}/log/bwa_cov_track_data.{fasta_prefix}.{pe_datatype}.{track_name}.fixmate.log",
        sort="{fasta_dir}/log/bwa_cov_track_data.{fasta_prefix}.{pe_datatype}.{track_name}.sort.log",
        markdup="{fasta_dir}/log/bwa_cov_track_data.{fasta_prefix}.{pe_datatype}.{track_name}.markdup.log",
        cluster_log="{fasta_dir}/log/bwa_cov_track_data.{fasta_prefix}.{pe_datatype}.{track_name}.cluster.log",
        cluster_err="{fasta_dir}/log/bwa_cov_track_data.{fasta_prefix}.{pe_datatype}.{track_name}.cluster.err"
    benchmark:
        "{fasta_dir}/log/bwa_cov_track_data.{fasta_prefix}.{pe_datatype}.{track_name}.benchmark.txt"

rule calculate_coverage:
    input:
        bam="{bam_dir}/{bam_prefix}.{datatype}.bam",
        csi="{bam_dir}/{bam_prefix}.{datatype}.bam.csi",
        log_dir=ancient("{bam_dir}/log/"),
    output:
        per_base="{bam_dir}/{bam_prefix}.{datatype}.{cov_settings}.per-base.bed.gz"
    params:
        min_mapq= lambda wildcards: parse_option("min_mapping_quality", parameters["tool_options"]["mosdepth"]["options"][wildcards.cov_settings][wildcards.datatype], " -Q ", none_value=0),
        blacklist_flags= lambda wildcards: parse_option("blacklist_flags", parameters["tool_options"]["mosdepth"]["options"][wildcards.cov_settings], " -F "),
    log:
        std="{bam_dir}/log/calculate_coverage.{bam_prefix}.{datatype}.{cov_settings}.log",
        cluster_log="{bam_dir}/log/calculate_coverage.{bam_prefix}.{datatype}.{cov_settings}.cluster.log",
        cluster_err="{bam_dir}/log/calculate_coverage.{bam_prefix}.{datatype}.{cov_settings}.cluster.err"
    benchmark:
        "{bam_dir}/log/calculate_coverage.{bam_prefix}.{cov_settings}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("calculate_coverage"),
        cpus=parameters["threads"]["mosdepth"],
        time=parameters["time"]["mosdepth"],
        mem=parameters["memory_mb"]["mosdepth"]
    threads: parameters["threads"]["mosdepth"]

    shell:
        " PREFIX={output.per_base};"
        " PREFIX=${{PREFIX%.per-base.bed.gz}}; "
        " mosdepth -t {threads} {params.blacklist_flags} {params.min_mapq} ${{PREFIX}} {input.bam} 2>{log.std}; "

rule create_coverage_table:
    input:
        per_base=rules.calculate_coverage.output.per_base,
        log_dir=ancient("{bam_dir}/log/"),
    output:
        stat_file="{bam_dir}/{bam_prefix}.{datatype}_{cov_settings}.win{window}.step{step}.stat",
        all_stat_file="{bam_dir}/{bam_prefix}.{datatype}_{cov_settings}.win{window}.step{step}.all.stat"
    log:
        std="{bam_dir}/log/create_coverage_table.{bam_prefix}.{datatype}.{cov_settings}.{window}.{step}.log",
        cluster_log="{bam_dir}/log/create_coverage_table.{bam_prefix}.{datatype}.{cov_settings}.{window}.{step}.cluster.log",
        cluster_err="{bam_dir}/log/create_coverage_table.{bam_prefix}.{datatype}.{cov_settings}.{window}.{step}.cluster.err"
    benchmark:
        "{bam_dir}/log/create_coverage_table.{bam_prefix}.{datatype}.{cov_settings}.{window}.{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_coverage_table"),
        cpus=parameters["threads"]["create_coverage_table"],
        time=parameters["time"]["create_coverage_table"],
        mem=parameters["memory_mb"]["create_coverage_table"]
    threads: parameters["threads"]["create_coverage_table"]

    shell:
        " PREFIX={output.stat_file};"
        " PREFIX=${{PREFIX%.stat}}; "
        " get_windows_stats_mosdepth_per_base_file.py -i {input.per_base} -w {wildcards.window} -s {wildcards.step} "
        "     -c bed -o ${{PREFIX}} 2>{log.std}; "

rule create_bedgraph_from_coverage_table:
    input:
        stat_file= rules.create_coverage_table.output.stat_file,
        log_dir=ancient("{bam_dir}/log/"),
    output:
        bedgraph="{bam_dir}/{bam_prefix}.{datatype}_{cov_settings}_{cov_type}_{track_type, coverage}.win{window}.step{step}.track.bedgraph"
    params:
        coverage_col= lambda wildcards: 6 if wildcards.cov_type == "mean" else 7
    log:
        std="{bam_dir}/log/create_bedgraph_from_coverage_table.{bam_prefix}.{datatype}.{cov_settings}.{cov_type}.{track_type}.{window}.{step}.log",
        cluster_log="{bam_dir}/log/create_bedgraph_from_coverage_table.{bam_prefix}.{datatype}.{cov_settings}.{cov_type}.{track_type}.{window}.{step}.cluster.log",
        cluster_err="{bam_dir}/log/create_bedgraph_from_coverage_table.{bam_prefix}.{datatype}.{cov_settings}.{cov_type}.{track_type}.{window}.{step}.cluster.err"
    benchmark:
        "{bam_dir}/log/create_bedgraph_from_coverage_table.{bam_prefix}.{datatype}.{cov_settings}.{cov_type}.{track_type}.{window}.{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_bedgraph_from_coverage_table"),
        cpus=parameters["threads"]["create_bedgraph_from_coverage_table"],
        time=parameters["time"]["create_bedgraph_from_coverage_table"],
        mem=parameters["memory_mb"]["create_bedgraph_from_coverage_table"]
    threads: parameters["threads"]["create_bedgraph_from_coverage_table"]

    shell:
        " tail -n +2 {input.stat_file} | cut -f 1,2,3,{params.coverage_col}  > {output.bedgraph} 2>{log.std}; "

rule draw_coverage_heatmap:
    input:
        stat_file="{fasta_dir}/{fasta_prefix}/alignment/{fasta_prefix}.{datatype}_{cov_settings}.win{window}.step{step}.stat",
        whitelist="{fasta_dir}/{fasta_prefix}.{scaffold_length}.whitelist",
        orderlist="{fasta_dir}/{fasta_prefix}.{scaffold_length}.orderlist",
        len_file="{fasta_dir}/{fasta_prefix}.len",
        all_stat_file="{fasta_dir}/{fasta_prefix}/alignment/{fasta_prefix}.{datatype}_{cov_settings}.win{window}.step{step}.all.stat",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        png="{fasta_dir}/assembly_qc/trackplots/{fasta_prefix}/{fasta_prefix}.{datatype}.{track_type, coverage}_{cov_settings}.{scaffold_length}.win{window}.step{step}.png",
        split_png="{fasta_dir}/assembly_qc/trackplots/{fasta_prefix}/{fasta_prefix}.{datatype}.{track_type, coverage}_{cov_settings}.{scaffold_length}.win{window}.step{step}.split_thresholds.png"

    log:
        std="{fasta_dir}/log/draw_coverage_heatmap.{fasta_prefix}.{datatype}.{track_type}.{cov_settings}.{scaffold_length}.{window}.{step}.log",
        cluster_log="{fasta_dir}/log/draw_coverage_heatmap.{fasta_prefix}.{track_type}.{datatype}.{cov_settings}.{scaffold_length}.{window}.{step}.cluster.log",
        cluster_err="{fasta_dir}/log/draw_coverage_heatmap.{fasta_prefix}.{track_type}.{datatype}.{cov_settings}.{scaffold_length}.{window}.{step}.cluster.err"
    benchmark:
        "{fasta_dir}/log/draw_coverage_heatmap.{fasta_prefix}.{datatype}.{track_type}.{cov_settings}.{scaffold_length}.{window}.{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("draw_coverage_heatmap"),
        cpus=parameters["threads"]["draw_coverage_heatmap"],
        time=parameters["time"]["draw_coverage_heatmap"],
        mem=parameters["memory_mb"]["draw_coverage_heatmap"]
    threads: parameters["threads"]["draw_coverage_heatmap"]

    shell:
        " PREFIX={output.png}; "
        " PREFIX=${{PREFIX%.png}}; "
        " draw_coverage.py -i {input.stat_file} -a {input.whitelist}  -w {wildcards.window} -s {wildcards.step} "
        "     -z {input.orderlist}  -n {input.len_file} -m `tail -n 1 {input.all_stat_file} | cut -f 6` --stranded_end  "
        "     --hide_track_label  --coverage_column_name_list median --rounded -o ${{PREFIX}} > {log.std} 2>&1; "
        " draw_coverage.py -i {input.stat_file} -a {input.whitelist}  -w {wildcards.window} -s {wildcards.step} "
        "     -z {input.orderlist}  -n {input.len_file} -m `tail -n 1 {input.all_stat_file} | cut -f 6` --stranded_end  "
        "     --hide_track_label  --split_coverage_thresholds --coverage_column_name_list median --rounded "
        "     -o ${{PREFIX}}.split_thresholds >> {log.std} 2>&1; "

rule gather_coverage_track_files:
    input:
        bedgraph="{fasta_dir}/{fasta_prefix}/alignment/{fasta_prefix}.{datatype}_{cov_settings}_{cov_type}_{track_type}.win{window}.step{step}.track.bedgraph",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.{datatype}_{cov_settings}_{cov_type}_{track_type, coverage}.win{window}.step{step}.track.bedgraph",
    log:
        std="{fasta_dir}/log/gather_coverage_track_files.{fasta_prefix}.{datatype}.{track_type}.{cov_settings}.{cov_type}.{window}.{step}.log",
        cluster_log="{fasta_dir}/log/gather_coverage_track_files.{fasta_prefix}.{track_type}.{datatype}.{cov_settings}.{cov_type}.{window}.{step}.cluster.log",
        cluster_err="{fasta_dir}/log/gather_coverage_track_files.{fasta_prefix}.{track_type}.{datatype}.{cov_settings}.{cov_type}.{window}.{step}.cluster.err"
    benchmark:
        "{fasta_dir}/log/gather_coverage_track_files.{fasta_prefix}.{datatype}.{track_type}.{cov_settings}.{cov_type}.{window}.{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_links"),
        cpus=parameters["threads"]["create_links"],
        time=parameters["time"]["create_links"],
        mem=parameters["memory_mb"]["create_links"]
    threads: parameters["threads"]["create_links"]

    shell:
        " cp -f {input.bedgraph} `dirname {output.bedgraph}` > {log.std} 2>&1; "

rule create_per_hap_coverage_tracks_for_merged_haplotype:
    input:
        bedgraph="{fasta_dir}/{genome_prefix}.{assembly_stage}.{haplotype}/alignment/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}_{cov_settings}_{cov_type}_{track_type}.win{window}.step{step}.track.bedgraph",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        bedgraph="{fasta_dir}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{merged_haplotype}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.{haplotype, hap[^/.@_]*}@{datatype}_{cov_settings}_{cov_type}_{track_type, coverage}.win{window}.step{step}.track.bedgraph",
    log:
        std="{fasta_dir}/log/create_per_hap_coverage_tracks_for_merged_haplotype.{genome_prefix}.{assembly_stage}.{merged_haplotype}.{haplotype}.{datatype}.{track_type}.{cov_settings}.{cov_type}.{window}.{step}.log",
        cluster_log="{fasta_dir}/log/create_per_hap_coverage_tracks_for_merged_haplotype.{genome_prefix}.{assembly_stage}.{merged_haplotype}.{haplotype}.{datatype}.{track_type}.{cov_settings}.{cov_type}.{window}.{step}.cluster.log",
        cluster_err="{fasta_dir}/log/create_per_hap_coverage_tracks_for_merged_haplotype.{genome_prefix}.{assembly_stage}.{merged_haplotype}.{haplotype}.{datatype}.{track_type}.{cov_settings}.{cov_type}.{window}.{step}.cluster.err"
    benchmark:
        "{fasta_dir}/log/gather_coverage_track_files.{genome_prefix}.{assembly_stage}.{merged_haplotype}.{haplotype}.{datatype}.{track_type}.{cov_settings}.{cov_type}.{window}.{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_per_hap_coverage_tracks_for_merged_haplotype"),
        cpus=parameters["threads"]["create_per_hap_coverage_tracks_for_merged_haplotype"],
        time=parameters["time"]["create_per_hap_coverage_tracks_for_merged_haplotype"],
        mem=parameters["memory_mb"]["create_per_hap_coverage_tracks_for_merged_haplotype"]
    threads: parameters["threads"]["create_per_hap_coverage_tracks_for_merged_haplotype"]

    shell:
        " sed 's/^/{wildcards.haplotype}./' {input.bedgraph} > {output.bedgraph} 2>{log.std}; "
