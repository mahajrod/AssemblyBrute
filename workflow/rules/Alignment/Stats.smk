
rule bam_stats:
    input:
        bam="{bam_prefix}.bam"
    output:
        stats="{bam_prefix}.bam.stats"
        #fai_alias=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta.fai"
    log:
        std="{bam_prefix}.bam.stats.log",
        cluster_log="{bam_prefix}.bam_stats.cluster.log",
        cluster_err="{bam_prefix}.bam_stats.cluster.err"
    benchmark:
        "{bam_prefix}.bam_stats.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["bam_stats"] ,
        time=parameters["time"]["bam_stats"],
        mem=parameters["memory_mb"]["bam_stats"]
    threads: parameters["threads"]["bam_stats"]

    shell:
        " samtools stats -@ {threads} --insert-size {}  {input.bam} > {output.stats} >{log.std} 2>&1;"

