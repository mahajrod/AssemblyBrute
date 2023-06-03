ruleorder: minimap2_cov > minimap2_purge_dups_reads

rule minimap2_cov: # TODO: add nanopore support
    input:
        fastq=lambda wildcards: expand(output_dict["data"] / ("fastq/%s/filtered/{fileprefix}%s" % (wildcards.datatype, config["fastq_extension"])),
                     fileprefix=input_file_prefix_dict[wildcards.datatype],
                     allow_missing=True),
        reference=out_dir_path  / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.fasta"
    output:
        bam=out_dir_path  / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.{datatype}.bam"
        #paf=out_dir_path  / ("purge_dups/{assembler}/{haplotype}/%s.purge_dups.{assembler}.{haplotype}.minimap2.{fileprefix}.paf.gz" % config["genome_name"])
    params:
        index_size=lambda wildcards: parse_option("index_size", parameters["tool_options"]["minimap2"][wildcards.datatype], " -I "),
        alignment_scheme=lambda wildcards: parse_option("alignment_scheme", parameters["tool_options"]["minimap2"][wildcards.datatype], " -x "),
        sort_threads=parameters["threads"]["samtools_sort"],
        minimap_threads=parameters["threads"]["minimap2"],
        per_thread_sort_mem=parameters["memory_mb"]["samtools_sort"],
    log:
        minimap2=output_dict["log"]  / "minimap2_cov.{prev_stage_parameters}.{curation_parameters}.{haplotype}.{genome_prefix}.{datatype}.minimap2.log",
        sort=output_dict["log"]  / "minimap2_cov.{prev_stage_parameters}.{curation_parameters}.{haplotype}.{genome_prefix}.{datatype}.sort.log",
        index=output_dict["log"]  / "minimap2_cov.{prev_stage_parameters}.{curation_parameters}.{haplotype}.{genome_prefix}.{datatype}.index.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_cov.{prev_stage_parameters}.{curation_parameters}.{haplotype}.{genome_prefix}.{datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_cov.{prev_stage_parameters}.{curation_parameters}.{haplotype}.{genome_prefix}.{datatype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_cov.{prev_stage_parameters}.{curation_parameters}.{haplotype}.{genome_prefix}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["minimap2"] + parameters["threads"]["samtools_sort"],
        time=parameters["time"]["minimap2"],
        mem=parameters["memory_mb"]["minimap2"] + (parameters["memory_mb"]["samtools_sort"] * parameters["threads"]["samtools_sort"])
    threads: parameters["threads"]["minimap2"] + parameters["threads"]["samtools_sort"]

    shell:
        " TMPDIR=`dirname {output.bam}`; "
        " minimap2 {params.alignment_scheme} {params.index_size} -a -t {params.minimap_threads}  {input.reference}  "
        " {input.fastq} 2>{log.minimap2} |  samtools sort -T ${{TMPDIR}} -@ {params.sort_threads} "
        " -m {params.per_thread_sort_mem}M -o {output.bam} 2>{log.sort};"
        " samtools index -@ {threads} {output.bam} > {log.index} 2>&1 "

rule calculate_coverage:
    input:
        bam=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.{datatype}.bam"
    output:
        per_base=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.{datatype}.per-base.bed.gz"
    log:
        std=output_dict["log"]  / "calculate_coverage.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{datatype}.log",
        cluster_log=output_dict["cluster_log"] / "calculate_coverage.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "calculate_coverage.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{datatype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "calculate_coverage.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["mosdepth"],
        time=parameters["time"]["mosdepth"],
        mem=parameters["memory_mb"]["mosdepth"]
    threads: parameters["threads"]["mosdepth"]

    shell:
        " PREFIX={output.per_base};"
        " PREFIX=${{PREFIX%.per-base.bed.gz}}; "
        " mosdepth -t {threads} ${{PREFIX}} {input.bam} 2>{log.std}; "

rule create_coverage_track:
    input:
        per_base=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.{datatype}.per-base.bed.gz"
    output:
        stat_file=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.{datatype}.stat"
    params:
        bin_size=lambda wildcards: stage_dict["curation"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.curation_parameters]["option_set"]["bin_size"],
        step_size=lambda wildcards: stage_dict["curation"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.curation_parameters]["option_set"]["bin_size"]
    log:
        std=output_dict["log"]  / "create_coverage_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{datatype}.log",
        cp=output_dict["log"]  / "create_coverage_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{datatype}.cp.log",
        cluster_log=output_dict["cluster_log"] / "create_coverage_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{datatype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_coverage_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{datatype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_coverage_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["create_coverage_track"],
        time=parameters["time"]["create_coverage_track"],
        mem=parameters["memory_mb"]["create_coverage_track"]
    threads: parameters["threads"]["create_coverage_track"]

    shell:
        " PREFIX={output.stat_file};"
        " PREFIX=${{PREFIX%.stat}}; "
        " get_windows_stats_mosdepth_per_base_file.py -i {input.per_base} -w {params.bin_size} -s {params.step_size} "
        " -o ${{PREFIX}} 2>{log.std}; "
        " cp ${{PREFIX}}.win{params.bin_size}.step{params.step_size}.stat ${{PREFIX}}.stat > {log.cp} 2>&1;"

