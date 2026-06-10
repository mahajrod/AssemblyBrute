ruleorder: get_track_stats > create_coverage_table
ruleorder: create_bedgraph_from_coverage_table > create_bedgraph_track
if "purge_dups" in config["stage_list"]:
    ruleorder: minimap2_cov > minimap2_purge_dups_reads

#wildcard_constraints:
#    longread_datatype="|".join(config["long_read_data"])

rule minimap2_cov: # TODO: add nanopore support
    input:
        fastq=lambda wildcards: expand(output_dict["data"] / ("%s/%s/%s/{fileprefix}%s" % (datatype_format_dict[wildcards.longread_datatype],
                                                                                                 wildcards.longread_datatype,
                                                                                                 "filtered" if wildcards.longread_datatype in config["filtered_data"] else "raw",
                                                                                                 config[datatype_format_dict[wildcards.longread_datatype] + "_extension"])),
                     fileprefix=input_file_prefix_dict[wildcards.longread_datatype] if datatype_format_dict[wildcards.longread_datatype] == "fastq" else input_fasta_file_prefix_dict[wildcards.longread_datatype],
                     allow_missing=True),
        reference=ancient(out_dir_path  / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta")
    output:
        bam=out_dir_path  / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/{track_type, coverage}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^/]+}/{genome_prefix}.{assembly_stage}.{haplotype}.{longread_datatype}.bam"
        #paf=out_dir_path  / ("purge_dups/{assembler}/{haplotype}/%s.purge_dups.{assembler}.{haplotype}.minimap2.{fileprefix}.paf.gz" % config["genome_name"])
    params:
        index_size=lambda wildcards: parse_option("index_size", parameters["tool_options"]["minimap2"][wildcards.longread_datatype], " -I "),
        alignment_scheme=lambda wildcards: parse_option("alignment_scheme", parameters["tool_options"]["minimap2"][wildcards.longread_datatype], " -x "),
        sort_threads=parameters["threads"]["samtools_sort"],
        minimap_threads=parameters["threads"]["minimap2"],
        per_thread_sort_mem=parameters["memory_mb"]["samtools_sort"],
    log:
        minimap2=output_dict["log"]  / "minimap2_cov.{assembly_stage}.{parameters}.{track_type}.{haplotype}.{genome_prefix}.{longread_datatype}.minimap2.log",
        sort=output_dict["log"]  / "minimap2_cov.{assembly_stage}.{parameters}.{track_type}.{haplotype}.{genome_prefix}.{longread_datatype}.sort.log",
        index=output_dict["log"]  / "minimap2_cov.{assembly_stage}.{parameters}.{track_type}.{haplotype}.{genome_prefix}.{longread_datatype}.index.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_cov.{assembly_stage}.{parameters}.{track_type}.{haplotype}.{genome_prefix}.{longread_datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_cov.{assembly_stage}.{parameters}.{track_type}.{haplotype}.{genome_prefix}.{longread_datatype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_cov.{assembly_stage}.{parameters}.{track_type}.{haplotype}.{genome_prefix}.{longread_datatype}.benchmark.txt"
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
        " {input.fastq} 2>{log.minimap2} |  samtools sort -T ${{TMPDIR}} -@ {params.sort_threads} "
        " -m {params.per_thread_sort_mem}M -o {output.bam} 2>{log.sort};"
        #" samtools index -@ {threads} {output.bam} > {log.index} 2>&1 "


rule bwa_cov:
    input:
        forward_fastqs=lambda wildcards: expand(output_dict["data"] / ("%s/%s/%s/{pairprefix}%s%s" % (datatype_format_dict[wildcards.pe_datatype],
                                                                                                      wildcards.pe_datatype,
                                                                                                      "filtered" if wildcards.pe_datatype in config["filtered_data"] else "raw",
                                                                                                      "_1" if wildcards.datatype in config["filtered_data"] else input_reverse_suffix_dict[wildcards.pe_datatype],
                                                                                                      datatype_extension_dict[wildcards.pe_datatype])),
                     pairprefix=input_pairprefix_dict[wildcards.datatype],
                     allow_missing=True),
        reverse_fastqs=lambda wildcards: expand(output_dict["data"] / ("%s/%s/%s/{pairprefix}%s%s" % (datatype_format_dict[wildcards.pe_datatype],
                                                                                                      wildcards.pe_datatype,
                                                                                                      "filtered" if wildcards.pe_datatype in config["filtered_data"] else "raw",
                                                                                                      "_2" if wildcards.pe_datatype in config["filtered_data"] else input_reverse_suffix_dict[wildcards.pe_datatype],
                                                                                                      datatype_extension_dict[wildcards.pe_datatype])),
                     pairprefix=input_pairprefix_dict[wildcards.datatype],
                     allow_missing=True),
        reference=out_dir_path  / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
        reference_index=out_dir_path  / ("{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta%s" % (".bwt" if config["bwa_tool"] == "bwa" else ".bwt.2bit.64")),
    output:
        bam=out_dir_path  / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/{track_type, coverage}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^/]+}/{genome_prefix}.{assembly_stage}.{haplotype}.{pe_datatype, illumina}.bam"

    params:
        bwa_tool=config["bwa_tool"],
        bwa_threads=parameters["threads"]["bwa_map"],
        sort_threads=parameters["threads"]["samtools_sort"],
        fixmate_threads=parameters["threads"]["samtools_fixmate"],
        markdup_threads=parameters["threads"]["samtools_markdup"],
        per_thread_sort_mem=parameters["memory_mb"]["samtools_sort_per_thread"],
        genome_prefix=config["genome_prefix"]
    log:
        bwa=output_dict["log"]  / "minimap2_cov.{assembly_stage}.{parameters}.{track_type}.{haplotype}.{genome_prefix}.{pe_datatype}.bwa.log",
        fixmate=output_dict["log"]  / "minimap2_cov.{assembly_stage}.{parameters}.{track_type}.{haplotype}.{genome_prefix}.{pe_datatype}.fixmate.log",
        sort=output_dict["log"]  / "minimap2_cov.{assembly_stage}.{parameters}.{track_type}.{haplotype}.{genome_prefix}.{pe_datatype}.sort.log",
        markdup=output_dict["log"]  / "minimap2_cov.{assembly_stage}.{parameters}.{track_type}.{haplotype}.{genome_prefix}.{pe_datatype}.markdup.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_cov.{assembly_stage}.{parameters}.{track_type}.{haplotype}.{genome_prefix}.{pe_datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_cov.{assembly_stage}.{parameters}.{track_type}.{haplotype}.{genome_prefix}.{pe_datatype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_cov.{assembly_stage}.{parameters}.{track_type}.{haplotype}.{genome_prefix}.{pe_datatype}.benchmark.txt"
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
        " {params.bwa_tool} mem  -t {params.bwa_threads} {input.reference} <(gunzip -c {input.forward_fastqs}) <(gunzip -c {input.reverse_fastqs}) "
        " -R  \'@RG\\tID:{params.genome_prefix}\\tPU:x\\tSM:{params.genome_prefix}\\tPL:Illumina\\tLB:x\' 2>{log.bwa} | "
        " samtools fixmate -@ {params.fixmate_threads} -m - -  2>{log.fixmate} | "
        " samtools sort -T {{TMP_PREFIX}} -@ {params.sort_threads} -m {params.per_thread_sort_mem}M 2>{log.sort} | "
        " samtools markdup -@ {params.markdup_threads} - {output.bam} 2>{log.markdup}"

rule calculate_coverage:
    input:
        bam=ancient(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/{track_type}/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.bam"),
        #bai=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/{track_type}/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.bam.bai"
        csi=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/{track_type}/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.bam.csi"

    output:
        per_base=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/{track_type, coverage}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^./]+}/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.{settings}.per-base.bed.gz"
    params:
        min_mapq= lambda wildcards: parse_option("min_mapping_quality", parameters["tool_options"]["mosdepth"]["options"][wildcards.settings][wildcards.datatype], " -Q ", none_value=0),
        blacklist_flags= lambda wildcards: parse_option("blacklist_flags", parameters["tool_options"]["mosdepth"]["options"][wildcards.settings], " -F "),
    log:
        std=output_dict["log"]  / "calculate_coverage.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{datatype}.{settings}.log",
        cluster_log=output_dict["cluster_log"] / "calculate_coverage.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{datatype}.{settings}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "calculate_coverage.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{datatype}.{settings}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "calculate_coverage.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{datatype}.{settings}.benchmark.txt"
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
    output:
        stat_file=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/{track_type, coverage}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^./]+}/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype, [^./]+}_{settings, [^./]+}.win{window, [0-9]+}.step{step, [0-9]+}.stat",
        all_stat_file=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/{track_type, coverage}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^./]+}/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype, [^._/]+}_{settings, [^./]+}.win{window, [0-9]+}.step{step, [0-9]+}.all.stat"
    log:
        std=output_dict["log"]  / "create_coverage_table.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{datatype}.{settings}.{window}.{step}.log",
        cluster_log=output_dict["cluster_log"] / "create_coverage_table.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{datatype}.{settings}.{window}.{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_coverage_table.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{datatype}.{settings}.{window}.{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_coverage_table.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{datatype}.{settings}.{window}.{step}.benchmark.txt"
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
        " -c bed -o ${{PREFIX}} 2>{log.std}; "

rule create_bedgraph_from_coverage_table:
    input:
        stat_file=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/{track_type}/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}_{settings}.win{window}.step{step}.stat"
    output:
        bedgraph=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/tracks/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^./]+}/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype, [^._/]+}_{settings, [^./]+}_{cov_type, [^./]+}_{track_type, coverage}.win{window, [0-9]+}.step{step, [0-9]+}.track.bedgraph"
    params:
        coverage_col= lambda wildcards: 6 if wildcards.cov_type == "mean" else 7
    log:
        std=output_dict["log"]  / "create_bedgraph_from_coverage_table.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{datatype}.{settings}.{cov_type}.{window}.{step}.log",
        cluster_log=output_dict["cluster_log"] / "create_bedgraph_from_coverage_table.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{datatype}.{settings}.{cov_type}.{window}.{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_bedgraph_from_coverage_table.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{datatype}.{settings}.{cov_type}.{window}.{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_bedgraph_from_coverage_table.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{datatype}.{settings}.{cov_type}.{window}.{step}.benchmark.txt"
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
        stat_file=rules.create_coverage_table.output.stat_file,
        whitelist=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.{scaffold_length}.whitelist",
        orderlist=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.{scaffold_length}.orderlist",
        len_file=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
        all_stat_file=rules.create_coverage_table.output.all_stat_file
    output:
        png=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/trackplots/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^./]+}/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype, [^./]+}.{track_type, coverage}_{settings, [^./]+}.{scaffold_length}.win{window, [0-9]+}.step{step, [0-9]+}.png",
        split_png=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/trackplots/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^./]+}/{genome_prefix}.{assembly_stage}.{haplotype}.{datatype, [^./]+}.{track_type, coverage}_{settings, [^./]+}.{scaffold_length}.win{window, [0-9]+}.step{step, [0-9]+}.split_thresholds.png"

    log:
        std=output_dict["log"]  / "draw_coverage_heatmap.{assembly_stage}..{parameters}..{track_type}.{scaffold_length}.{genome_prefix}.{haplotype}.{datatype}.{settings}.{window}.{step}.log",
        cluster_log=output_dict["cluster_log"] / "draw_coverage_heatmap.{assembly_stage}..{parameters}..{track_type}.{scaffold_length}.{genome_prefix}.{haplotype}.{datatype}.{settings}.{window}.{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "draw_coverage_heatmap.{assembly_stage}..{parameters}..{track_type}.{scaffold_length}.{genome_prefix}.{haplotype}.{datatype}.{settings}.{window}.{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "draw_coverage_heatmap.{assembly_stage}..{parameters}..{track_type}.{scaffold_length}.{genome_prefix}.{haplotype}.{datatype}.{settings}.{window}.{step}.benchmark.txt"
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
        " -z {input.orderlist}  -n {input.len_file} -m `tail -n 1 {input.all_stat_file} | cut -f 6` --stranded_end  "
        " --hide_track_label  --coverage_column_name_list median --rounded -o ${{PREFIX}} > {log.std} 2>&1; "
        " draw_coverage.py -i {input.stat_file} -a {input.whitelist}  -w {wildcards.window} -s {wildcards.step} "
        " -z {input.orderlist}  -n {input.len_file} -m `tail -n 1 {input.all_stat_file} | cut -f 6` --stranded_end  "
        " --hide_track_label  --split_coverage_thresholds --coverage_column_name_list median --rounded "
        "-o ${{PREFIX}}.split_thresholds >> {log.std} 2>&1; "