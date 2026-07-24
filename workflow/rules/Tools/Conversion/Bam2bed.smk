
rule bam2bed:
    input:
        bam="{bam_dir}/{bam_prefix}.rmdup.bam",
        log_dir=ancient("{bam_dir}/log/")
    output:
        bed="{bam_dir}/{bam_prefix}.rmdup.bed"
    log:
        convert="{bam_dir}/log/bam2bed.{bam_prefix}.convert.log",
        sort="{bam_dir}/log/bam2bed.{bam_prefix}.sort.log",
        cluster_log="{bam_dir}/log/bam2bed.{bam_prefix}.cluster.log",
        cluster_err="{bam_dir}/log/bam2bed.{bam_prefix}.cluster.err"
    benchmark:
        "{bam_dir}/log/bam2bed.{bam_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("bam_to_bed"),
        cpus=parameters["threads"]["bam2bed"] ,
        time=parameters["time"]["bam2bed"],
        mem=parameters["memory_mb"]["bam2bed"]
    threads: parameters["threads"]["bam2bed"]

    shell:
        " bamToBed -i {input.bam} 2>{log.convert} | sort -S{resources.mem}M -k 4 > {output.bed} 2>{log.sort}"
