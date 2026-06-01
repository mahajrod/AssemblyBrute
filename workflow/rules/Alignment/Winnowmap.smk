

rule winnowmap: #
    input:
        fastq=lambda wildcards: expand(output_dict["data"] / ("%s/%s/%s/{fileprefix}%s" % (datatype_format_dict[wildcards.datatype],
                                                                                           wildcards.datatype,
                                                                                           "filtered" if wildcards.datatype in config["filtered_data"] else "raw",
                                                                                           config[datatype_format_dict[wildcards.datatype] + "_extension"])),
                     fileprefix=input_file_prefix_dict[wildcards.datatype] if datatype_format_dict[wildcards.datatype] == "fastq" else input_fasta_file_prefix_dict[wildcards.datatype],
                     allow_missing=True),
        reference="{directory}/{fasta_prefix}.fasta",
        repetitive_kmers=directory("{directory}/kmer/meryl/{fasta_prefix}.{kmer_length}.meryl.repetitive"),
        log_dir="{directory}/log/",
        benchmark_dir="{directory}/benchmark/"
    output:
        #bam=out_dir_path  / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{fileprefix}.bwa.bam"
        bam="{directory}/alignment/{datatype}/winnowmap/{fasta_prefix}.k{kmer_length}.{datatype}.winnowmap.bam"
    params:
        sort_threads=parameters["threads"]["samtools_sort"],
        minimap_threads=parameters["threads"]["minimap2"],
        per_thread_sort_mem=parameters["memory_mb"]["samtools_sort"],
    log:
        map="{directory}/log/winnowmap.{datatype}.{fasta_prefix}.{kmer_length}.{datatype}.map.log",
        sort="{directory}/log/winnowmap.{datatype}.{fasta_prefix}.{kmer_length}.{datatype}.sort.log",
        cluster_log="{directory}/log/winnowmap.{datatype}.{fasta_prefix}.{kmer_length}.{datatype}.cluster.log",
        cluster_err="{directory}/log/winnowmap.{datatype}.{fasta_prefix}.{kmer_length}.{datatype}.cluster.err"
    benchmark:
        "{directory}/benchmark/winnowmap.{datatype}.{fasta_prefix}.{kmer_length}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("minimap2_cov"),
        cpus=parameters["threads"]["minimap2"] + parameters["threads"]["samtools_sort"],
        time=parameters["time"]["minimap2"],
        mem=parameters["memory_mb"]["minimap2"] + (parameters["memory_mb"]["samtools_sort"] * parameters["threads"]["samtools_sort"])
    threads: parameters["threads"]["minimap2"] + parameters["threads"]["samtools_sort"]
    shell:
        " winnowmap -t {threads} -W {input.repetitive_kmers} -ax map-pb {input.reference} {input.fastq} 2>{log.map} | "
        " samtools sort -@ {params.sort_threads}  -m {params.per_thread_sort_mem}M -o {output.bam} - 2>{log.sort}; "