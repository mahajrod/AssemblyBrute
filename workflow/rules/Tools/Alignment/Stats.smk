
rule bam_stats:
    input:
        bam="{bam_dir}/{bam_prefix}.bam",
        log_dir=ancient("{bam_dir}/log/")
    output:
        stats="{bam_dir}/{bam_prefix}.bam.general_stats"
    params:
        max_insert_size=parameters["tool_options"]["samtools_stats"]["max_insert_size"]
    log:
        std="{bam_dir}/log/bam.stats.{bam_prefix}.log",
        cluster_log="{bam_dir}/log/bam_stats.cluster.{bam_prefix}.log",
        cluster_err="{bam_dir}/log/bam_stats.cluster.{bam_prefix}.err"
    benchmark:
        "{bam_dir}/log/bam_stats.{bam_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("bam_stats"),
        cpus=parameters["threads"]["bam_stats"] ,
        time=parameters["time"]["bam_stats"],
        mem=parameters["memory_mb"]["bam_stats"],
        bam_stats=1
    threads: parameters["threads"]["bam_stats"]

    shell:
        " samtools stats -@ {threads} --insert-size {params.max_insert_size}  {input.bam} > {output.stats} 2>{log.std}; "
