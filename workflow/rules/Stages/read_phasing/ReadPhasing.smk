ruleorder: extract_pe_reads_by_unique_hap_kmers > extract_se_reads_by_unique_hap_kmers

rule meryl_extract_unique_hap_kmers:
    input:
        target_hap_db_dir=config["out_dir"]  / "{stage}/{parameters}/{genome_prefix}.{stage}.{haplotype}/kmer/meryl/{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}.meryl",
        rest_hap_db_dirs=lambda wildcards: expand(config["out_dir"]  / ("%s/%s/%s.%s.{haplotype}/kmer/meryl/%s.%s.{haplotype}.%s.meryl" % (wildcards.stage,
                                                                                                                                           wildcards.parameters,
                                                                                                                                           wildcards.genome_prefix,
                                                                                                                                           wildcards.stage,
                                                                                                                                           wildcards.genome_prefix,
                                                                                                                                           wildcards.stage,
                                                                                                                                           wildcards.assembly_kmer_length)) ,
                                                 haplotype=set(stage_dict[wildcards.stage].parameters[wildcards.parameters]["haplotype_list"]) - {wildcards.haplotype},
                                                 allow_missing=True)
    output:
        unique_hap_db_dir=directory(config["out_dir"]  / "{stage}/{parameters}/{genome_prefix}.{stage}.{haplotype}/kmer/meryl/{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length, [^.]+}.unique.meryl"),
    log:
        std=config["out_dir"]/ "log/meryl_extract_unique_hap_kmers.{genome_prefix}.{stage}.{haplotype}.{parameters}.{assembly_kmer_length}.log",
        cluster_log=config["out_dir"] / "log/meryl_extract_unique_hap_kmers.{genome_prefix}.{stage}.{haplotype}.{parameters}.{assembly_kmer_length}.cluster.log",
        cluster_err=config["out_dir"] / "log/meryl_extract_unique_hap_kmers.{genome_prefix}.{stage}.{haplotype}.{parameters}.{assembly_kmer_length}.cluster.err"
    benchmark:
        config["out_dir"]/ "log/meryl_extract_unique_hap_kmers.{genome_prefix}.{stage}.{haplotype}.{parameters}.{assembly_kmer_length}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("meryl_extract_unique_hap_kmers"),
        cpus=parameters["threads"]["meryl_extract_unique_hap_kmers"],
        time=parameters["time"]["meryl_extract_unique_hap_kmers"],
        mem=parameters["memory_mb"]["meryl_extract_unique_hap_kmers"],
    threads:
        parameters["threads"]["meryl_extract_unique_hap_kmers"]
    shell:
         " meryl threads={threads} memory={resources.mem}m difference {input.target_hap_db_dir} {input.rest_hap_db_dirs} output {output.unique_hap_db_dir} > {log.std} 2>&1; "

use rule merge_meryl as merge_alien_kmers with:
    input:
        lambda wildcards: expand(config["out_dir"]  / ("%s/%s/%s.%s.{haplotype}/kmer/meryl/%s.%s.{haplotype}.%s.unique.meryl" % (wildcards.stage,
                                                                                                                                 wildcards.parameters,
                                                                                                                                 wildcards.genome_prefix,
                                                                                                                                 wildcards.stage,
                                                                                                                                 wildcards.genome_prefix,
                                                                                                                                 wildcards.stage,
                                                                                                                                 wildcards.assembly_kmer_length)),
                                                 haplotype=set(stage_dict[wildcards.stage].parameters[wildcards.parameters]["haplotype_list"]) - {wildcards.haplotype},
                                                 allow_missing=True),

    output:
        db_dir=directory(config["out_dir"] / "{stage}/{parameters}/{genome_prefix}.{stage}.{haplotype}/kmer/meryl/{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}.alien.unique.meryl"),
    log:
        count_log=config["out_dir"] / "log/merge_alien_kmers.{parameters}.{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}.count.log",
        cluster_log=config["out_dir"] / "log/merge_alien_kmers.{parameters}.{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}.cluster.log",
        cluster_err=config["out_dir"] / "log/merge_alien_kmers.{parameters}.{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}.cluster.err"
    benchmark:
        config["out_dir"] / "log/merge_alien_kmers.{parameters}.{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}.benchmark.txt"

rule extract_pe_reads_by_unique_hap_kmers:
    input:
        alien_kmer_db_dir=config["out_dir"] / "{stage}/{parameters}/{genome_prefix}.{stage}.{haplotype}/kmer/meryl/{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}.alien.unique.meryl",
        forward_read=lambda wildcards: config["out_dir"]  / ("data/{0}/{1}/{2}{3}{4}".format(wildcards.datatype,
                                                                                                "final",
                                                                                                wildcards.pairprefix,
                                                                                                config["data"][wildcards.datatype]["conv_fwd_sfx"],
                                                                                                config["data"][wildcards.datatype]["conv_ext"] )),

        reverse_read=lambda wildcards: config["out_dir"] / ("data/{0}/{1}/{2}{3}{4}".format(wildcards.datatype,
                                                                                                "final",
                                                                                                wildcards.pairprefix,
                                                                                                config["data"][wildcards.datatype]["conv_rev_sfx"],
                                                                                                config["data"][wildcards.datatype]["conv_ext"] )),
    output:
        forward_hap_read=config["out_dir"] / ("{stage}/{parameters}/{genome_prefix}.{stage}.{haplotype}/reads/{datatype}/{assembly_kmer_length}/{pairprefix}%s.fastq.gz" % config["fwd_fastq_sfx"]),
        reverse_hap_read=config["out_dir"] / ("{stage}/{parameters}/{genome_prefix}.{stage}.{haplotype}/reads/{datatype}/{assembly_kmer_length}/{pairprefix}%s.fastq.gz" % config["rev_fastq_sfx"]),
    log:
        std=config["out_dir"] / "log/extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{pairprefix}.{assembly_kmer_length}.{genome_prefix}.{haplotype}.log",
        cluster_log=config["out_dir"] / "log/extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{pairprefix}.{assembly_kmer_length}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{pairprefix}.{assembly_kmer_length}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{pairprefix}.{assembly_kmer_length}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("extract_pe_reads_by_unique_hap_kmers"),
        cpus=parameters["threads"]["extract_reads_by_unique_hap_kmers"],
        time=parameters["time"]["extract_reads_by_unique_hap_kmers"],
        mem=parameters["memory_mb"]["extract_reads_by_unique_hap_kmers"],
    threads:
        parameters["threads"]["extract_reads_by_unique_hap_kmers"]
    shell:
         " meryl-lookup -exclude -sequence {input.forward_read} {input.reverse_read} "
         "     -mers {input.alien_kmer_db_dir} -output {output.forward_hap_read} {output.reverse_hap_read} > {log.std} 2>&1;"

rule extract_se_reads_by_unique_hap_kmers:
    input:
        alien_kmer_db_dir=config["out_dir"] / "{stage}/{parameters}/{genome_prefix}.{stage}.{haplotype}/kmer/meryl/{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}.alien.unique.meryl",
        se_read=lambda wildcards: config["out_dir"]  / ("data/{0}/{1}/{2}{3}".format(wildcards.datatype,
                                                                                     "final",
                                                                                     wildcards.fileprefix,
                                                                                     config["data"][wildcards.datatype]["conv_ext"])),
    output:
        hap_se_read=config["out_dir"] / "{stage}/{parameters}/{genome_prefix}.{stage}.{haplotype}/reads/{datatype}/{assembly_kmer_length}/{fileprefix}.fastq.gz",
    log:
        std=config["out_dir"] / "log/extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{genome_prefix}.{haplotype}.log",
        cluster_log=config["out_dir"] / "log/extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("extract_se_reads_by_unique_hap_kmers"),
        cpus=parameters["threads"]["extract_reads_by_unique_hap_kmers"],
        time=parameters["time"]["extract_reads_by_unique_hap_kmers"],
        mem=parameters["memory_mb"]["extract_reads_by_unique_hap_kmers"],
    threads:
        parameters["threads"]["extract_reads_by_unique_hap_kmers"]
    shell:
         " meryl-lookup -exclude -sequence {input.se_read} "
         "     -mers {input.alien_kmer_db_dir} -output {output.hap_se_read} > {log.std} 2>&1;"

use rule extract_se_reads_by_unique_hap_kmers as extract_se_fasta_reads_by_unique_hap_kmers with:
    input:
        alien_kmer_db_dir=config["out_dir"] / "{stage}/{parameters}/{genome_prefix}.{stage}.{haplotype}/kmer/meryl/{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}.alien.unique.meryl",
        se_read=lambda wildcards: config["out_dir"]  / ("data/{0}/{1}/{2}{3}".format(wildcards.datatype,
                                                                                     "final",
                                                                                     wildcards.fileprefix,
                                                                                     config["data"][wildcards.datatype]["conv_ext"])),
    output:
        hap_se_read=config["out_dir"] / "{stage}/{parameters}/{genome_prefix}.{stage}.{haplotype}/reads/{datatype}/{assembly_kmer_length}/{fileprefix}.fasta.gz",


"""
rule extract_se_reads_from_fasta_by_unique_hap_kmers: #TODO: merge with extract_se_reads_by_unique_hap_kmers:
    input:
        rest_hap_db_dirs=lambda wildcards: expand(config["out_dir"]  / ("%s/%s/%s.%s.{haplotype}/kmer/meryl/%s.%s.{haplotype}.%s.unique.meryl" % (wildcards.stage,
                                                                                                              wildcards.parameters,
                                                                                                              config["genome_prefix"],
                                                                                                              wildcards.stage,
                                                                                                              config["genome_prefix"],
                                                                                                              wildcards.stage,
                                                                                                              wildcards.assembly_kmer_length)),
                                                 haplotype=set(stage_dict[wildcards.stage].parameters[wildcards.parameters]["haplotype_list"]) - {wildcards.haplotype},
                                                 allow_missing=True),
        se_read=lambda wildcards: config["out_dir"]  / ("data/fasta/{0}/{1}/{2}{3}".format(wildcards.datatype,
                                                                                        "filtered" if wildcards.datatype in config["filtered_data"] else "raw",
                                                                                        wildcards.fileprefix,
                                                                                        config["fasta_extension"])),
    output:
        hap_se_read=config["out_dir"] / ("{stage}/{parameters}/%s.{stage}.{haplotype}/reads/{datatype}/{assembly_kmer_length}/{fileprefix}.fasta.gz" % config["genome_prefix"]),
    log:
        std=config["out_dir"] / "log/extract_se_reads_from_fasta_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{haplotype}.log",
        cluster_log=config["out_dir"] / "log/extract_se_reads_from_fasta_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{haplotype}.cluster.log",
        cluster_err=config["out_dir"] / "log/extract_se_reads_from_fasta_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/extract_se_reads_from_fasta_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("extract_se_reads_from_fasta_by_unique_hap_kmers"),
        cpus=parameters["threads"]["extract_reads_by_unique_hap_kmers"],
        time=parameters["time"]["extract_reads_by_unique_hap_kmers"],
        mem=parameters["memory_mb"]["extract_reads_by_unique_hap_kmers"],
    threads:
        parameters["threads"]["extract_reads_by_unique_hap_kmers"]
    shell:
         " meryl-lookup -exclude -sequence {input.se_read} "
         " -mers {input.rest_hap_db_dirs} -output {output.hap_se_read} > {log.std} 2>&1; "
"""