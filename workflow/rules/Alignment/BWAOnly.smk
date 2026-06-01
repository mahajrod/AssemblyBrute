localrules: link_bwa_only_bam

rule bwa_map: #
    input:
        index=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta.ann",
        reference=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
        #fastq=lambda wildcards: output_dict["data"] / "fastq/hic/raw/{0}{1}".format(wildcards.fileprefix, config["fastq_extension"]) if wildcards.phasing_kmer_length == "NA" else \
        #                        out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}".format(config["phasing_stage"], #wildcards.assembly_stage,
        #                                                                                 detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
        #                                                                                 wildcards.haplotype,
        #                                                                                 wildcards.phasing_kmer_length,
        #                                                                                 wildcards.fileprefix,
        #                                                                                 config["fastq_extension"]
        #                                                                                 ),
        forward_fastq=lambda wildcards: output_dict["data"] / "fastq/hic/{0}/{1}{2}{3}".format("filtered" if "hic" in config["filtered_data"] else "raw",
                                                                                               wildcards.pairprefix,
                                                                                               "_1" if "hic" in config["filtered_data"] else input_forward_suffix_dict["hic"],
                                                                                               config["fastq_extension"]) if wildcards.phasing_kmer_length == "NA" else \
                                out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}{6}".format(config["phasing_stage"], #wildcards.assembly_stage,
                                                                                            detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
                                                                                            wildcards.haplotype,
                                                                                            wildcards.phasing_kmer_length,
                                                                                            wildcards.pairprefix,
                                                                                            "_1",
                                                                                            config["fastq_extension"]),
        reverse_fastq=lambda wildcards: output_dict["data"] / "fastq/hic/{0}/{1}{2}{3}".format("filtered" if "hic" in config["filtered_data"] else "raw",
                                                                                               wildcards.pairprefix,
                                                                                               "_2" if "hic" in config["filtered_data"] else input_reverse_suffix_dict["hic"],
                                                                                               config["fastq_extension"]) if wildcards.phasing_kmer_length == "NA" else \
                                out_dir_path / "{0}/{1}/fastq/{2}/{3}/hic/{4}{5}{6}".format(config["phasing_stage"], #wildcards.assembly_stage,
                                                                                            detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."), #wildcards.parameters,
                                                                                            wildcards.haplotype,
                                                                                            wildcards.phasing_kmer_length,
                                                                                            wildcards.pairprefix,
                                                                                            "_2",
                                                                                            config["fastq_extension"]),
    output:
        #bam=out_dir_path  / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{fileprefix}.bwa.bam"
        raw_bam=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{pairprefix, [^/]+}.bwa.raw.bam",
        bam=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{pairprefix, [^/]+}.bwa.bam"
    params:
        id="{0}_hic".format(config["genome_prefix"]),
        bwa_tool=config["bwa_tool"]
    log:
        map=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.map.log",
        sort=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.sort.log",
        #filter=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.filter.log",
        cluster_log=output_dict["cluster_log"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("bwa_map"),
        cpus=parameters["threads"]["bwa_map"] ,
        time=parameters["time"]["bwa_map"],
        mem=parameters["memory_mb"]["bwa_map"]
    threads: parameters["threads"]["bwa_map"]
    shell:
        " {params.bwa_tool} mem -SP5M -t {threads} -R  \'@RG\\tID:{params.id}\\tPU:x\\tSM:{params.id}\\tPL:illumina\\tLB:x\' "
        " {input.reference} {input.forward_fastq} {input.reverse_fastq} 2>{log.map} | samtools view -Sb - > {output.bam} 2>{log.sort} "

rule link_bwa_only_bam: #
    input:
        raw_bam=out_dir_path / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{pairprefix}.bwa.raw.bam",
    output:
        bam=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/{haplotype, [^.]+}/alignment/{phasing_kmer_length, [^/]+}/{genome_prefix, [^/]+}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{pairprefix, [^/]+}.bwa.bam"

    log:
        log=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.log",
        cluster_log=output_dict["cluster_log"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_fastq_links"),
        cpus=parameters["threads"]["create_fastq_links"] ,
        time=parameters["time"]["create_fastq_links"],
        mem=parameters["memory_mb"]["create_fastq_links"]
    threads: parameters["threads"]["create_fastq_links"]
    shell:
        " ln -sf `basename {input.raw_bam}` {output.bam} 2>{log.log} "

