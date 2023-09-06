rule bwa_index:
    input:
        fasta="{fasta_prefix}.fasta"
    output:
        index="{fasta_prefix}.fasta%s" % (".bwt" if config["bwa_tool"] == "bwa" else ".bwt.2bit.64"), #  or (config["other_tool_option_sets"]["mapping_pipeline"] == "arima")
        index_ann="{fasta_prefix}.fasta.ann"
    params:
        bwa_tool=config["bwa_tool"] # if config["other_tool_option_sets"]["mapping_pipeline"] != "arima" else "bwa",
    log:
        std="{fasta_prefix}.bwa_index.log",
        cluster_log="{fasta_prefix}.bwa_index.cluster.log",
        cluster_err="{fasta_prefix}.bwa_index.cluster.err"
    benchmark:
        "{fasta_prefix}.bwa_index.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("bwa_index"),
        cpus=parameters["threads"]["bwa_index"] ,
        time=parameters["time"]["bwa_index"],
        mem=parameters["memory_mb"]["bwa_index"]
    threads: parameters["threads"]["bwa_index"]

    shell:
        " {params.bwa_tool} index {input.fasta} 1>{log.std} 2>&1;"

"""
rule bwa_index:
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta"
    output:
        index=out_dir_path / ("{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta%s" % (".bwt" if config["bwa_tool"] == "bwa" else ".bwt.2bit.64")), #  or (config["other_tool_option_sets"]["mapping_pipeline"] == "arima")
    params:
        bwa_tool=config["bwa_tool"] # if config["other_tool_option_sets"]["mapping_pipeline"] != "arima" else "bwa",
    log:
        std=output_dict["log"]  / "bwa_index.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "bwa_index.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_index.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_indexs.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        cpus=parameters["threads"]["bwa_index"] ,
        time=parameters["time"]["bwa_index"],
        mem=parameters["memory_mb"]["bwa_index"]
    threads: parameters["threads"]["bwa_index"]

    shell:
        " {params.bwa_tool} index {input.fasta} 1>{log.std} 2>&1;"
"""

rule ref_faidx:
    input:
        fasta="{fasta_prefix}.fasta"
    output:
        fai="{fasta_prefix}.fasta.fai",
        #fai_alias=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta.fai"
    log:
        std="{fasta_prefix}.ref_faidx.log",
        cluster_log="{fasta_prefix}.ref_faidx.cluster.log",
        cluster_err="{fasta_prefix}.ref_faidx.cluster.err"
    benchmark:
        "{fasta_prefix}.ref_faidx.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("ref_faidx"),
        cpus=parameters["threads"]["ref_faidx"] ,
        time=parameters["time"]["ref_faidx"],
        mem=parameters["memory_mb"]["ref_faidx"]
    threads: parameters["threads"]["ref_faidx"]

    shell:
        " samtools faidx -o {output.fai} {input.fasta} > {log.std} 2>&1;"

rule ref_dict:
    input:
        fasta="{fasta_prefix}.fasta"
    output:
        dict="{fasta_prefix}.dict"
    log:
        std="{fasta_prefix}.ref_dict.log",
        cluster_log="{fasta_prefix}.ref_dict.cluster.log",
        cluster_err="{fasta_prefix}.ref_dict.cluster.err"
    benchmark:
        "{fasta_prefix}.ref_dict.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("ref_dict"),
        cpus=parameters["threads"]["ref_dict"] ,
        time=parameters["time"]["ref_dict"],
        mem=parameters["memory_mb"]["ref_dict"]
    threads: parameters["threads"]["ref_dict"]

    shell:
        " picard CreateSequenceDictionary R={input.fasta} O={output.dict} > {log.std} 2>&1;"

rule index_bam:
    input:
        bam="{bam_prefix}.bam"
    output:
        bai="{bam_prefix}.bam.bai"
    log:
        std="{bam_prefix}.index_bam.log",
        cluster_log="{bam_prefix}.index_bam.cluster.log",
        cluster_err="{bam_prefix}.index_bam.cluster.err"
    benchmark:
        "{bam_prefix}.index_bam.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("index_bam"),
        cpus=parameters["threads"]["samtools_index"] ,
        time=parameters["time"]["samtools_index"],
        mem=parameters["memory_mb"]["samtools_index"]
    threads: parameters["threads"]["samtools_index"]
    shell:
        " samtools index -@ {threads} {input} > {log.std} 2>&1; "
