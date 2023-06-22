rule meryl_assembly:
    input:
        fasta=out_dir_path / "{stage}/{parameters}/{genome_prefix}.{stage}.{haplotype}.fasta"
    output:
        db_dir=directory(out_dir_path / "{stage}/{parameters}/kmer/{genome_prefix}.{stage}.{haplotype, [^./]+}.{assembly_kmer_length, [^./]+}")
    log:
        std=output_dict["log"] / "meryl_assembly.{genome_prefix}.{stage}.{parameters}.{haplotype}.{assembly_kmer_length}.log",
        cluster_log=output_dict["cluster_log"] / "meryl_assembly.{genome_prefix}.{stage}.{parameters}.{haplotype}.{assembly_kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "meryl_assembly.{genome_prefix}.{stage}.{parameters}.{haplotype}.{assembly_kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "meryl_assembly.{genome_prefix}.{stage}.{parameters}.{haplotype}.{assembly_kmer_length}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["meryl_assembly"],
        time=parameters["time"]["meryl_assembly"],
        mem=parameters["memory_mb"]["meryl_assembly"],
    threads:
        parameters["threads"]["meryl_assembly"]
    shell:
         " meryl k={wildcards.assembly_kmer_length} threads={threads} memory={resources.mem}m count "
         " output {output.db_dir} {input} 1>{log.std} 2>&1;"


rule meryl_extract_unique_hap_kmers:
    input:
        target_hap_db_dir=out_dir_path / "{stage}/{parameters}/kmer/{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}",
        rest_hap_db_dirs=lambda wildcards: expand(out_dir_path / ("%s/%s/kmer/%s.%s.{haplotype}.%s" % (wildcards.stage,
                                                                                                       wildcards.parameters,
                                                                                                       wildcards.genome_prefix,
                                                                                                       wildcards.stage,
                                                                                                       wildcards.assembly_kmer_length)) ,
                                                 haplotype=set(stage_dict[wildcards.stage]["parameters"][wildcards.parameters]["haplotype_list"]) - {wildcards.haplotype},
                                                 allow_missing=True)
    output:
        unique_hap_db_dir=directory(out_dir_path / "{stage}/{parameters}/kmer/{genome_prefix}.{stage}.{haplotype, [^./]+}.{assembly_kmer_length, [^.]+}.unique"),
        #unique_hap2_db_dir=out_dir_path / "{stage}/{parameters}/kmer/{genome_prefix}.{stage}.hap2.{assembly_kmer_length}.unique"
    log:
        std=output_dict["log"] / "meryl_extract_unique_hap_kmers.{genome_prefix}.{stage}.{parameters}.{assembly_kmer_length}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "meryl_extract_unique_hap_kmers.{genome_prefix}.{stage}.{parameters}.{assembly_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "meryl_extract_unique_hap_kmers.{genome_prefix}.{stage}.{parameters}.{assembly_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "meryl_extract_unique_hap_kmers.{genome_prefix}.{stage}.{parameters}.{assembly_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["meryl_extract_unique_hap_kmers"],
        time=parameters["time"]["meryl_extract_unique_hap_kmers"],
        mem=parameters["memory_mb"]["meryl_extract_unique_hap_kmers"],
    threads:
        parameters["threads"]["meryl_extract_unique_hap_kmers"]
    shell:
         " meryl threads={threads} memory={resources.mem}m difference {input.target_hap_db_dir} {input.rest_hap_db_dirs} output {output.unique_hap_db_dir} > {log.std} 2>&1; "

rule extract_pe_reads_by_unique_hap_kmers:
    input:
        rest_hap_db_dirs=lambda wildcards: expand(out_dir_path / ("%s/%s/kmer/%s.%s.{haplotype}.%s.unique" % (wildcards.stage,
                                                                                                              wildcards.parameters,
                                                                                                              config["genome_prefix"],
                                                                                                              wildcards.stage,
                                                                                                              wildcards.assembly_kmer_length)),
                                                 haplotype=set(stage_dict[wildcards.stage]["parameters"][wildcards.parameters]["haplotype_list"]) - {wildcards.haplotype},
                                                 allow_missing=True),
        forward_read=lambda wildcards: output_dict["data"]  / ("fastq/{0}/{1}/{2}{3}{4}".format(wildcards.datatype,
                                                                                                "filtered" if wildcards.datatype in config["filtered_data"] else "raw",
                                                                                                wildcards.pairprefix,
                                                                                                input_forward_suffix_dict[wildcards.datatype],
                                                                                                config["fastq_extension"])),

        reverse_read=lambda wildcards: output_dict["data"]  / ("fastq/{0}/{1}/{2}{3}{4}".format(wildcards.datatype,
                                                                                                "filtered" if wildcards.datatype in config["filtered_data"] else "raw",
                                                                                                wildcards.pairprefix,
                                                                                                input_reverse_suffix_dict[wildcards.datatype],
                                                                                                config["fastq_extension"])),
    output:
        forward_hap_read=out_dir_path / "{stage}/{parameters}/fastq/{haplotype, [^.]+}/{assembly_kmer_length, [^./]+}/{datatype, [^/]+}/{pairprefix, [^/]+}_1.fastq.gz", # TODO: change to forward_suffix
        reverse_hap_read=out_dir_path / "{stage}/{parameters}/fastq/{haplotype, [^.]+}/{assembly_kmer_length, [^./]+}/{datatype, [^/]+}/{pairprefix, [^/]+}_2.fastq.gz", # TODO: change to reverse_suffix
    log:
        std=output_dict["log"] / "extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{pairprefix}.{assembly_kmer_length}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{pairprefix}.{assembly_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{pairprefix}.{assembly_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{pairprefix}.{assembly_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["extract_reads_by_unique_hap_kmers"],
        time=parameters["time"]["extract_reads_by_unique_hap_kmers"],
        mem=parameters["memory_mb"]["extract_reads_by_unique_hap_kmers"],
    threads:
        parameters["threads"]["extract_reads_by_unique_hap_kmers"]
    shell:
         " meryl-lookup -exclude -sequence {input.forward_read} {input.reverse_read} "
         " -mers {input.rest_hap_db_dirs} -output {output.forward_hap_read} {output.reverse_hap_read} > {log.std} 2>&1;"

rule extract_se_reads_by_unique_hap_kmers:
    input:
        rest_hap_db_dirs=lambda wildcards: expand(out_dir_path / ("%s/%s/kmer/%s.%s.{haplotype}.%s.unique" % (wildcards.stage,
                                                                                                              wildcards.parameters,
                                                                                                              config["genome_prefix"],
                                                                                                              wildcards.stage,
                                                                                                              wildcards.assembly_kmer_length)),
                                                 haplotype=set(stage_dict[wildcards.stage]["parameters"][wildcards.parameters]["haplotype_list"]) - {wildcards.haplotype},
                                                 allow_missing=True),
        se_read=lambda wildcards: output_dict["data"]  / ("fastq/{0}/{1}/{2}{3}".format(wildcards.datatype,
                                                                                        "filtered" if wildcards.datatype in config["filtered_data"] else "raw",
                                                                                        wildcards.fileprefix,
                                                                                        config["fastq_extension"])),
    output:
        hap_se_read=out_dir_path / "{stage}/{parameters}/fastq/{haplotype, [^./]+}/{assembly_kmer_length, [^./]+}/{datatype, [^/]+}/{fileprefix, [^/]+}.fastq.gz",
    log:
        std=output_dict["log"] / "extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{haplotype}.log",
        #hap2=output_dict["log"] / "extract_reads_by_unique_hap_kmers.{stage}.{parameters}.{fileprefix}.{genome_prefix}.AK{assembly_kmer_length}.hap2.log",
        cluster_log=output_dict["cluster_log"] / "extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "extract_reads_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["extract_reads_by_unique_hap_kmers"],
        time=parameters["time"]["extract_reads_by_unique_hap_kmers"],
        mem=parameters["memory_mb"]["extract_reads_by_unique_hap_kmers"],
    threads:
        parameters["threads"]["extract_reads_by_unique_hap_kmers"]
    shell:
         " meryl-lookup -exclude -sequence {input.se_read} "
         " -mers {input.rest_hap_db_dirs} -output {output.hap_se_read} > {log.std} 2>&1;"

rule extract_se_reads_from_fasta_by_unique_hap_kmers: #TODO: merge with extract_se_reads_by_unique_hap_kmers:
    input:
        rest_hap_db_dirs=lambda wildcards: expand(out_dir_path / ("%s/%s/kmer/%s.%s.{haplotype}.%s.unique" % (wildcards.stage,
                                                                                                              wildcards.parameters,
                                                                                                              config["genome_prefix"],
                                                                                                              wildcards.stage,
                                                                                                              wildcards.assembly_kmer_length)),
                                                 haplotype=set(stage_dict[wildcards.stage]["parameters"][wildcards.parameters]["haplotype_list"]) - {wildcards.haplotype},
                                                 allow_missing=True),
        se_read=lambda wildcards: output_dict["data"]  / ("fasta/{0}/{1}/{2}{3}".format(wildcards.datatype,
                                                                                        "filtered" if wildcards.datatype in config["filtered_data"] else "raw",
                                                                                        wildcards.fileprefix,
                                                                                        config["fasta_extension"])),
    output:
        hap_se_read=out_dir_path / "{stage}/{parameters}/fasta/{haplotype, [^./]+}/{assembly_kmer_length, [^./]+}/{datatype, [^/]+}/{fileprefix, [^/]+}.fasta.gz",
    log:
        std=output_dict["log"] / "extract_se_reads_from_fasta_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{haplotype}.log",
        #hap2=output_dict["log"] / "extract_reads_by_unique_hap_kmers.{stage}.{parameters}.{fileprefix}.{genome_prefix}.AK{assembly_kmer_length}.hap2.log",
        cluster_log=output_dict["cluster_log"] / "extract_se_reads_from_fasta_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "extract_se_reads_from_fasta_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "extract_se_reads_from_fasta_by_unique_hap_kmers.{datatype}.{stage}.{parameters}.{fileprefix}.{assembly_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["extract_reads_by_unique_hap_kmers"],
        time=parameters["time"]["extract_reads_by_unique_hap_kmers"],
        mem=parameters["memory_mb"]["extract_reads_by_unique_hap_kmers"],
    threads:
        parameters["threads"]["extract_reads_by_unique_hap_kmers"]
    shell:
         " meryl-lookup -exclude -sequence {input.se_read} "
         " -mers {input.rest_hap_db_dirs} -output {output.hap_se_read} > {log.std} 2>&1;"