
rule bwa_index:
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta"
    output:
        index=out_dir_path / ("{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta%s" % (".bwt" if config["bwa_tool"] == "bwa" else ".bwt.2bit.64")),
    params:
        bwa_tool=config["bwa_tool"],
    log:
        std=output_dict["log"]  / "bwa_index.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "bwa_index.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_index.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_indexs.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["bwa_index"] ,
        time=parameters["time"]["bwa_index"],
        mem=parameters["memory_mb"]["bwa_index"]
    threads: parameters["threads"]["bwa_index"]

    shell:
        " {params.bwa_tool} index {input.fasta} 1>{log.std} 2>&1;"

rule ref_faidx:
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta"
    output:
        fai=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta.fai",
        #fai_alias=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta.fai"
    log:
        std=output_dict["log"]  / "ref_faidx.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "ref_faidx.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "ref_faidx.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "ref_faidx.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["ref_faidx"] ,
        time=parameters["time"]["ref_faidx"],
        mem=parameters["memory_mb"]["ref_faidx"]
    threads: parameters["threads"]["ref_faidx"]

    shell:
        " samtools faidx -o {output.fai} {input.fasta} >{log.std} 2>&1;"

rule ref_dict:
    input:
        fasta=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta"
    output:
        dict=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.dict"
    log:
        std=output_dict["log"]  / "ref_dict.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "ref_dict.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "ref_dict.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "ref_dict.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}..benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
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
        cpus=parameters["threads"]["samtools_index"] ,
        time=parameters["time"]["samtools_index"],
        mem=parameters["memory_mb"]["samtools_index"]
    threads: parameters["threads"]["samtools_index"]
    shell:
        " samtools index -@ {threads} {input} > {log.std} 2>&1; "
