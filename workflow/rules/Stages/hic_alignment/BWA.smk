
rule bwa_map: #
    input:
        index=config["out_dir"]  / "hic_alignment/{parameters}/{genome_prefix}.hic_alignment.{haplotype}.fasta.ann",
        reference=config["out_dir"]  / "hic_alignment/{parameters}/{genome_prefix}.hic_alignment.{haplotype}.fasta",
        forward_fastq=lambda wildcards: config["out_dir"] / "data//hic/final/{0}{1}{2}".format(wildcards.pairprefix,
                                                                                               config["data"]["hic"]["conv_fwd_sfx"],
                                                                                               config["data"]["hic"]["conv_ext"]) if wildcards.phasing_kmer_length == "NA" else \
                                config["out_dir"]  / "{0}/{1}/{2}.{0}.{3}/reads/hic/{4}/{5}{6}{7}".format(config["phasing_stage"],
                                                                                                    detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."),
                                                                                                    wildcards.genome_prefix,
                                                                                                    wildcards.haplotype,
                                                                                                    wildcards.phasing_kmer_length,
                                                                                                    wildcards.pairprefix,
                                                                                                    config["data"]["hic"]["conv_fwd_sfx"],
                                                                                                    config["data"]["hic"]["conv_ext"]),
        reverse_fastq=lambda wildcards: config["out_dir"] / "data//hic/final/{0}{1}{2}".format(wildcards.pairprefix,
                                                                                               config["data"]["hic"]["conv_rev_sfx"],
                                                                                               config["data"]["hic"]["conv_ext"]) if wildcards.phasing_kmer_length == "NA" else \
                                config["out_dir"]  / "{0}/{1}/{2}.{0}.{3}/reads/hic/{4}/{5}{6}{7}".format(config["phasing_stage"],
                                                                                                    detect_phasing_parameters(wildcards.parameters, config["phasing_stage"], stage_separator=".."),
                                                                                                    wildcards.genome_prefix,
                                                                                                    wildcards.haplotype,
                                                                                                    wildcards.phasing_kmer_length,
                                                                                                    wildcards.pairprefix,
                                                                                                    config["data"]["hic"]["conv_rev_sfx"],
                                                                                                    config["data"]["hic"]["conv_ext"]),
    output:
        bam=config["out_dir"] / "hic_alignment/{parameters, .*bwa_only.*|.*pairtools.*}/{genome_prefix}.hic_alignment.{haplotype, hap[^./]+}/alignment/{phasing_kmer_length}/{genome_prefix}.hic_alignment.{haplotype}.{phasing_kmer_length}.{pairprefix}.bwa.bam"
    log:
        map=config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.map.log",
        sort=config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.sort.log",
        cluster_log=config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/bwa_map.hic_alignment.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("bwa_map"),
        cpus=get_threads(parameters["threads"]["bwa_map"], "cpu"),
        time=parameters["time"]["bwa_map"],
        mem=parameters["memory_mb"]["bwa_map"]
    threads: parameters["threads"]["bwa_map"]
    shell:
        " bwa-mem2 mem -SP5M -t {threads} -R  \'@RG\\tID:{wildcards.genome_prefix}_hic\\tPU:x\\tSM:{wildcards.genome_prefix}_hic\\tPL:illumina\\tLB:x\' "
        "     {input.reference} {input.forward_fastq} {input.reverse_fastq} 2>{log.map} | samtools view -Sb - > {output.bam} 2>{log.sort} "
