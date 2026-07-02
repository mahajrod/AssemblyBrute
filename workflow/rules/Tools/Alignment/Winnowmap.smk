

rule winnowmap: #
    input:
        fastq=lambda wildcards: expand(config["out_dir"] / ("data/%s/final/{fileprefix}%s" % (wildcards.datatype,
                                                                                                    config["data"][wildcards.datatype]["conv_ext"])),
                     fileprefix=config["data"][wildcards.datatype]["conv_file_prefix_list"],
                     allow_missing=True),
        reference="{fasta_dir}/{fasta_prefix}.fasta",
        repetitive_kmers="{fasta_dir}/{fasta_prefix}/kmer/meryl/{fasta_prefix}.{kmer_length}.meryl.repetitive",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        bam="{fasta_dir}/alignment/{datatype}/winnowmap/{fasta_prefix}.k{kmer_length}.{datatype}.winnowmap.bam"
    params:
        sort_threads=parameters["threads"]["samtools_sort"],
        minimap_threads=parameters["threads"]["minimap2"],
        per_thread_sort_mem=parameters["memory_mb"]["samtools_sort"],
    log:
        map="{fasta_dir}/log/winnowmap.{datatype}.{fasta_prefix}.{kmer_length}.{datatype}.map.log",
        sort="{fasta_dir}/log/winnowmap.{datatype}.{fasta_prefix}.{kmer_length}.{datatype}.sort.log",
        cluster_log="{fasta_dir}/log/winnowmap.{datatype}.{fasta_prefix}.{kmer_length}.{datatype}.cluster.log",
        cluster_err="{fasta_dir}/log/winnowmap.{datatype}.{fasta_prefix}.{kmer_length}.{datatype}.cluster.err"
    benchmark:
        "{fasta_dir}/log/winnowmap.{datatype}.{fasta_prefix}.{kmer_length}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("minimap2_cov"),
        cpus=get_threads(parameters["threads"]["winnowmap"] + parameters["threads"]["samtools_sort"], "cpu"),
        time=parameters["time"]["winnowmap"],
        mem=parameters["memory_mb"]["winnowmap"] + (parameters["memory_mb"]["samtools_sort"] * parameters["threads"]["samtools_sort"])
    threads: parameters["threads"]["winnowmap"] + parameters["threads"]["samtools_sort"]
    shell:
        " winnowmap -t {threads} -W {input.repetitive_kmers} -ax map-pb {input.reference} {input.fastq} 2>{log.map} | "
        "     samtools sort -@ {params.sort_threads}  -m {params.per_thread_sort_mem}M -o {output.bam} - 2>{log.sort}; "