
rule rmdup:
    input:
        bam="{bam_dir}/{bam_prefix}.bwa.bam",
        stats="{bam_dir}/{bam_prefix}.bwa.bam.general_stats",
        log_dir="{bam_dir}/log/"
    output:
        bam="{bam_dir}/{bam_prefix, ^(?!.*rmdup).*}.rmdup.bam",
    params:
        sort_threads=parameters["threads"]["samtools_sort"],
        collate_threads=parameters["threads"]["samtools_collate"],
        fixmate_threads=parameters["threads"]["samtools_fixmate"],
        markdup_threads=parameters["threads"]["samtools_markdup"],
        sort_per_thread=parameters["memory_mb"]["samtools_sort_per_thread"]
    log:
        collate="{bam_dir}/log/rmdup.{bam_prefix}.collate.log",
        fixmate="{bam_dir}/log/rmdup.{bam_prefix}.fixmate.log",
        sort="{bam_dir}/log/rmdup.{bam_prefix}.sort.log",
        markdup="{bam_dir}/log/rmdup.{bam_prefix}.markdup.log",
        cluster_log="{bam_dir}/log/rmdup.{bam_prefix}.cluster.log",
        cluster_err="{bam_dir}/log/rmdup.{bam_prefix}.cluster.err"
    benchmark:
        "{bam_dir}/log/rmdup.{bam_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("rmdup"),
        cpus=parameters["threads"]["samtools_sort"] + parameters["threads"]["samtools_collate"] + parameters["threads"]["samtools_fixmate"] + parameters["threads"]["samtools_markdup"],
        time=parameters["time"]["rmdup"],
        mem=50000 + parameters["memory_mb"]["samtools_collate"] + parameters["memory_mb"]["samtools_fixmate"] + parameters["memory_mb"]["samtools_markdup"] + parameters["memory_mb"]["samtools_sort_per_thread"] * parameters["threads"]["samtools_sort"]
    threads: parameters["threads"]["samtools_sort"] + parameters["threads"]["samtools_collate"] + parameters["threads"]["samtools_fixmate"] + parameters["threads"]["samtools_markdup"]
    shell:
        " TMP_DIR=`dirname {output.bam}`; "
        " samtools collate -T ${{TMP_DIR}}/tmp.collate  -@ {params.collate_threads}  -O {input.bam} 2>{log.collate} | "
        " samtools fixmate -@ {params.fixmate_threads} -m - -  2>{log.fixmate} | "
        " samtools sort -T ${{TMP_DIR}}/tmp.sort -@ {params.sort_threads} -m {params.sort_per_thread}M 2>{log.sort} | "
        " samtools markdup -@ {params.markdup_threads} - {output.bam} > {log.markdup} 2>&1; "
