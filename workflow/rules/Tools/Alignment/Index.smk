
rule bwa_mem2_index:
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/")
    output:
        index_bwt="{fasta_dir}/{fasta_prefix}.fasta.bwt.2bit.64",
        index_ann="{fasta_dir}/{fasta_prefix}.fasta.ann",
        index_amb="{fasta_dir}/{fasta_prefix}.fasta.amb",
        index_0123="{fasta_dir}/{fasta_prefix}.fasta.0123",
        index_pac="{fasta_dir}/{fasta_prefix}.fasta.pac",
    log:
        std="{fasta_dir}/log/bwa_index.{fasta_prefix}.log",
        cluster_log="{fasta_dir}/log/bwa_index.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/log/bwa_index.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/log/{fasta_prefix}.bwa_index.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("bwa_index"),
        cpus=parameters["threads"]["bwa_index"] ,
        time=parameters["time"]["bwa_index"],
        mem=partial(get_memory, start_mem=parameters["memory_mb"]["bwa_index"], coeff=1.5, mode="exp")
    threads: parameters["threads"]["bwa_index"]

    shell:
        " bwa-mem2 index {input.fasta} 1>{log.std} 2>&1;"

rule ref_faidx:
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/")
    output:
        fai="{fasta_dir}/{fasta_prefix}.fasta.fai",
    log:
        std="{fasta_dir}/log/ref_faidx.{fasta_prefix}.log",
        cluster_log="{fasta_dir}/log/ref_faidx.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/log/ref_faidx.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/log/ref_faidx.{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("ref_faidx"),
        cpus=parameters["threads"]["ref_faidx"] ,
        time=parameters["time"]["ref_faidx"],
        mem=partial(get_memory, start_mem=parameters["memory_mb"]["ref_faidx"], coeff=1.5, mode="exp")
    threads: parameters["threads"]["ref_faidx"]

    shell:
        " samtools faidx -o {output.fai} {input.fasta} > {log.std} 2>&1;"

rule ref_dict:
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/")
    output:
        dict="{fasta_dir}/{fasta_prefix}.dict"
    log:
        std="{fasta_dir}/log/ref_dict.{fasta_prefix}.log",
        cluster_log="{fasta_dir}/log/ref_dict.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/log/ref_dict.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/log/{fasta_prefix}.ref_dict.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("ref_dict"),
        cpus=parameters["threads"]["ref_dict"] ,
        time=parameters["time"]["ref_dict"],
        mem=partial(get_memory, start_mem=parameters["memory_mb"]["ref_dict"], coeff=1.5, mode="exp")
    threads: parameters["threads"]["ref_dict"]

    shell:
        " picard CreateSequenceDictionary R={input.fasta} O={output.dict} > {log.std} 2>&1;"

rule index_bam:
    input:
        bam="{bam_dir}/{bam_prefix}.bam",
        log_dir=ancient("{bam_dir}/log/")
    output:
        bai="{bam_dir}/{bam_prefix}.bam.bai"
    log:
        std="{bam_dir}/log/index_bam.{bam_prefix}.log",
        cluster_log="{bam_dir}/log/index_bam.{bam_prefix}.cluster.log",
        cluster_err="{bam_dir}/log/index_bam.{bam_prefix}.cluster.err"
    benchmark:
        "{bam_dir}/log/index_bam.{bam_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("index_bam"),
        cpus=parameters["threads"]["samtools_index"] ,
        time=parameters["time"]["samtools_index"],
        mem=partial(get_memory, start_mem=parameters["memory_mb"]["samtools_index"], coeff=1.5, mode="exp")
    threads: parameters["threads"]["samtools_index"]
    shell:
        " samtools index -@ {threads} {input.bam} > {log.std} 2>&1; "

rule index_bam_csi:
    input:
        bam="{bam_dir}/{bam_prefix}.bam",
        log_dir=ancient("{bam_dir}/log/")
    output:
        bai="{bam_dir}/{bam_prefix}.bam.csi"
    log:
        std="{bam_dir}/log/index_bam.{bam_prefix}.log",
        cluster_log="{bam_dir}/log/index_bam.cluster.{bam_prefix}.log",
        cluster_err="{bam_dir}/log/index_bam.cluster.{bam_prefix}.err"
    benchmark:
        "{bam_dir}/log/index_bam.{bam_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("index_bam"),
        cpus=parameters["threads"]["samtools_index"] ,
        time=parameters["time"]["samtools_index"],
        mem=partial(get_memory, start_mem=parameters["memory_mb"]["samtools_index"], coeff=1.5, mode="exp")
    threads: parameters["threads"]["samtools_index"]
    shell:
        " samtools index -c -@ {threads} {input.bam} > {log.std} 2>&1; "
