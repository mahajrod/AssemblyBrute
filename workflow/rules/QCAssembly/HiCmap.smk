rule bwa_map_for_hic_map: #
    input:
        index="{fasta_dir}/{fasta_prefix}.fasta.ann",
        reference="{fasta_dir}/{fasta_prefix}.fasta",
        forward_fastq=lambda wildcards: ancient(config["out_dir"] / "data/hic/final/{0}{1}{2}".format(wildcards.pairprefix,
                                                                                                          config["data"]["hic"]["conv_fwd_sfx"],
                                                                                                          config["data"]["hic"]["conv_ext"])),
        reverse_fastq=lambda wildcards: ancient(config["out_dir"] / "data/hic/final/{0}{1}{2}".format(wildcards.pairprefix,
                                                                                                          config["data"]["hic"]["conv_rev_sfx"],
                                                                                                          config["data"]["hic"]["conv_ext"])),

        log_dir=ancient("{fasta_dir}/log/")
    output:
        bam=temp("{fasta_dir}/{fasta_prefix, .*combined.*|.*reordered.*}/alignment/NA/{fasta_prefix}.NA.{pairprefix}.bwa.bam")
    params:
        id="{0}_hic".format(config["genome_prefix"]),
    log:
        map="{fasta_dir}/log/bwa_map.{fasta_prefix}.NA.{pairprefix}.map.log",
        sort="{fasta_dir}/log/bwa_map.{fasta_prefix}.NA.{pairprefix}.sort.log",
        cluster_log="{fasta_dir}/log/bwa_map.{fasta_prefix}.NA.{pairprefix}.cluster.log",
        cluster_err="{fasta_dir}/log/bwa_map.{fasta_prefix}.NA.{pairprefix}.cluster.err"
    benchmark:
        "{fasta_dir}/log/bwa_map.{fasta_prefix}.NA.{pairprefix}.benchmark.txt"
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
        " bwa-mem2 mem  -T 0 -SP5 -t {threads} -R  \'@RG\\tID:{params.id}\\tPU:x\\tSM:{params.id}\\tPL:illumina\\tLB:x\' "
        "     {input.reference} {input.forward_fastq} {input.reverse_fastq} 2>{log.map} | samtools view -Sb - > {output.bam} 2>{log.sort} "