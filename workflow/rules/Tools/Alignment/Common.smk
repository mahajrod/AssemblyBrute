
rule rmdup:
    input:
        bam="{bam_dir}/{bam_prefix}.bwa.bam",
        stats="{bam_dir}/{bam_prefix}.bwa.bam.general_stats",
        log_dir=ancient("{bam_dir}/log/")
    output:
        bam="{bam_dir}/{bam_prefix}.rmdup.bam",
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
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("rmdup"),
        cpus=parameters["threads"]["samtools_sort"] + parameters["threads"]["samtools_collate"] + parameters["threads"]["samtools_fixmate"] + parameters["threads"]["samtools_markdup"],
        time=parameters["time"]["rmdup"],
        mem=50000 + parameters["memory_mb"]["samtools_collate"] + parameters["memory_mb"]["samtools_fixmate"] + parameters["memory_mb"]["samtools_markdup"] + parameters["memory_mb"]["samtools_sort_per_thread"] * parameters["threads"]["samtools_sort"]
    threads: parameters["threads"]["samtools_sort"] + parameters["threads"]["samtools_collate"] + parameters["threads"]["samtools_fixmate"] + parameters["threads"]["samtools_markdup"]
    shell:
        " TMP_DIR=`dirname {output.bam}`; "
        " samtools collate -T ${{TMP_DIR}}/tmp.collate  -@ {params.collate_threads}  -O {input.bam} 2>{log.collate} | "
        "     samtools fixmate -@ {params.fixmate_threads} -m - -  2>{log.fixmate} | "
        "     samtools sort -T ${{TMP_DIR}}/tmp.sort -@ {params.sort_threads} -m {params.sort_per_thread}M 2>{log.sort} | "
        "     samtools markdup -@ {params.markdup_threads} - {output.bam} > {log.markdup} 2>&1; "


if "hic" in config["data"]:
    rule bam_merge_hic_files:
        input:
            bams=expand("{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.{pairprefix}.bwa.bam",
                        allow_missing=True,
                        pairprefix=config["data"]["hic"]["pair_prefix_list"]),
            reference_fai="{fasta_dir}/{fasta_prefix}.fasta.fai",
            reference="{fasta_dir}/{fasta_prefix}.fasta",
            log_dir=ancient("{fasta_dir}/log/"),
        output:
            bam=temp("{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.bwa.bam")
        params:
            sort_threads=parameters["threads"]["samtools_sort"]
        log:
            std="{fasta_dir}/log/bam_merge_files.{fasta_prefix}.{phasing_kmer_length}.log",
            cluster_log="{fasta_dir}/log/bam_merge_files.{fasta_prefix}.{phasing_kmer_length}.cluster.log",
            cluster_err="{fasta_dir}/log/bam_merge_files.{fasta_prefix}.{phasing_kmer_length}.cluster.err"
        benchmark:
            "{fasta_dir}/log/bam_merge_files.{fasta_prefix}.{phasing_kmer_length}.benchmark.txt"
        conda:
            config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
        resources:
            queue=config["queue"]["cpu"]["name"],
            node_options=parse_node_list("bwa_merge_files"),
            cpus=parameters["threads"]["samtools_sort"] ,
            time=parameters["time"]["samtools_sort"],
            mem=parameters["memory_mb"]["samtools_sort"]
        threads: parameters["threads"]["samtools_sort"]
        shell:
            " samtools merge -@ {params.sort_threads} --no-PG -o {output.bam} {input.bams} >{log.std} 2>&1"