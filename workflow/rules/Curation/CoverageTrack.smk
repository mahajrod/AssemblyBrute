ruleorder: get_track_stats > create_coverage_table
ruleorder: create_bedgraph_from_coverage_table > create_bedgraph_track
if "purge_dups" in config["stage_list"]:
    ruleorder: minimap2_cov > minimap2_purge_dups_reads

rule minimap2_cov: # TODO: add nanopore support
    input:
        fastq=lambda wildcards: expand(output_dict["data"] / ("%s/%s/%s/{fileprefix}%s" % (datatype_format_dict[wildcards.datatype],
                                                                                           wildcards.datatype,
                                                                                           "filtered" if wildcards.datatype in config["filtered_data"] else "raw",
                                                                                           config[datatype_format_dict[wildcards.datatype] + "_extension"])),
                     fileprefix=input_file_prefix_dict[wildcards.datatype] if datatype_format_dict[wildcards.datatype] == "fastq" else input_fasta_file_prefix_dict[wildcards.datatype],
                     allow_missing=True),
        reference=out_dir_path  / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.fasta"
    output:
        bam=out_dir_path  / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype}.bam"
        #paf=out_dir_path  / ("purge_dups/{assembler}/{haplotype}/%s.purge_dups.{assembler}.{haplotype}.minimap2.{fileprefix}.paf.gz" % config["genome_name"])
    params:
        index_size=lambda wildcards: parse_option("index_size", parameters["tool_options"]["minimap2"][wildcards.datatype], " -I "),
        alignment_scheme=lambda wildcards: parse_option("alignment_scheme", parameters["tool_options"]["minimap2"][wildcards.datatype], " -x "),
        sort_threads=parameters["threads"]["samtools_sort"],
        minimap_threads=parameters["threads"]["minimap2"],
        per_thread_sort_mem=parameters["memory_mb"]["samtools_sort"],
    log:
        minimap2=output_dict["log"]  / "minimap2_cov.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.minimap2.log",
        sort=output_dict["log"]  / "minimap2_cov.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.sort.log",
        index=output_dict["log"]  / "minimap2_cov.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.index.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_cov.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_cov.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_cov.{prev_stage_parameters}.{curation_parameters}.{seq_type}.{haplotype}.{genome_prefix}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("minimap2_cov"),
        cpus=parameters["threads"]["minimap2"] + parameters["threads"]["samtools_sort"],
        time=parameters["time"]["minimap2"],
        mem=parameters["memory_mb"]["minimap2"] + (parameters["memory_mb"]["samtools_sort"] * parameters["threads"]["samtools_sort"])
    threads: parameters["threads"]["minimap2"] + parameters["threads"]["samtools_sort"]

    shell:
        " TMPDIR=`dirname {output.bam}`; "
        " minimap2 {params.alignment_scheme} {params.index_size} -a -t {params.minimap_threads}  {input.reference}  "
        " {input.fastq} 2>{log.minimap2} |  samtools sort -T ${{TMPDIR}} -@ {params.sort_threads} "
        " -m {params.per_thread_sort_mem}M -o {output.bam} 2>{log.sort};"
        #" samtools index -@ {threads} {output.bam} > {log.index} 2>&1 "

rule calculate_coverage:
    input:
        bam=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype}.bam",
        bai=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype}.bam.bai"
    output:
        per_base=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype}.per-base.bed.gz"
    params:
        min_mapq= lambda wildcards: parse_option("min_mapping_quality", parameters["tool_options"]["mosdepth"][wildcards.datatype], " -Q ", none_value=0),
    log:
        std=output_dict["log"]  / "calculate_coverage.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.log",
        cluster_log=output_dict["cluster_log"] / "calculate_coverage.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "calculate_coverage.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "calculate_coverage.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("calculate_coverage"),
        cpus=parameters["threads"]["mosdepth"],
        time=parameters["time"]["mosdepth"],
        mem=parameters["memory_mb"]["mosdepth"]
    threads: parameters["threads"]["mosdepth"]

    shell:
        " PREFIX={output.per_base};"
        " PREFIX=${{PREFIX%.per-base.bed.gz}}; "
        " mosdepth -t {threads} ${{PREFIX}} {input.bam} 2>{log.std}; "

rule create_coverage_table:
    input:
        per_base=rules.calculate_coverage.output.per_base,
    output:
        stat_file=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype}.win{window}.step{step}.stat",
        all_stat_file=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype}.win{window}.step{step}.all.stat"
    #params:
    #    bin_size=lambda wildcards: stage_dict["curation"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.curation_parameters]["option_set"]["bin_size"],
    #    step_size=lambda wildcards: stage_dict["curation"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.curation_parameters]["option_set"]["bin_size"]
    log:
        std=output_dict["log"]  / "create_coverage_table.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.{window}.{step}.log",
        #cp=output_dict["log"]  / "create_coverage_table.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{datatype}.{window}.{step}.cp.log",
        cluster_log=output_dict["cluster_log"] / "create_coverage_table.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.{window}.{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_coverage_table.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.{window}.{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_coverage_table.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.{window}.{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
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
        #" cp ${{PREFIX}}.win{params.bin_size}.step{params.step_size}.stat ${{PREFIX}}.stat > {log.cp} 2>&1;"

rule create_bedgraph_from_coverage_table:
    input:
        stat_file=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype}.win{window}.step{step}.stat"
    output:
        bedgraph=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/{seq_type, [^./]+}/{genome_prefix}.input.{haplotype}.{datatype, [^./]+}_{cov_type, [^./]+}_coverage.win{window, [^./]+}.step{step, [^./]+}.track.bedgraph"
    params:
        coverage_col= lambda wildcards: 6 if wildcards.cov_type == "mean" else 7
    log:
        std=output_dict["log"]  / "create_bedgraph_from_coverage_table.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.{cov_type}.{window}.{step}.log",
        cluster_log=output_dict["cluster_log"] / "create_bedgraph_from_coverage_table.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.{cov_type}.{window}.{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_bedgraph_from_coverage_table.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.{cov_type}.{window}.{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_bedgraph_from_coverage_table.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.{cov_type}.{window}.{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("create_bedgraph_from_coverage_table"),
        cpus=parameters["threads"]["create_coverage_table"],
        time=parameters["time"]["create_coverage_table"],
        mem=parameters["memory_mb"]["create_coverage_table"]
    threads: parameters["threads"]["create_coverage_table"]

    shell:
        " tail -n +2 {input.stat_file} | cut -f 1,2,3,{params.coverage_col}  > {output.bedgraph} 2>{log.std}; "

rule draw_coverage_heatmap:
    input:
        stat_file=rules.create_coverage_table.output.stat_file,
        whitelist=rules.select_long_scaffolds.output.whitelist,
        orderlist=rules.select_long_scaffolds.output.orderlist,
        len_file=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.len",
        all_stat_file=rules.create_coverage_table.output.all_stat_file
    output:
        png=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^./]+}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype,  [^./]+}.coverage.win{window}.step{step}.png",
        split_png=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^./]+}/{seq_type}/{genome_prefix}.input.{haplotype}.{datatype, [^./]+}.coverage.win{window}.step{step}.split_thresholds.png"
    #params:
    #    bin_size=lambda wildcards: stage_dict["curation"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.curation_parameters]["option_set"]["bin_size"],
    #    step_size=lambda wildcards: stage_dict["curation"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.curation_parameters]["option_set"]["bin_size"]
    log:
        std=output_dict["log"]  / "draw_coverage_heatmap.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.{window}.{step}.log",
        cluster_log=output_dict["cluster_log"] / "draw_coverage_heatmap.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.{window}.{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "draw_coverage_heatmap.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.{window}.{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "draw_coverage_heatmap.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.{datatype}.{window}.{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
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