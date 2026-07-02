

rule yak_pe_reads:
    input:
        forward_fastqs=lambda wildcards: expand(output_dict["data"] / ("fastq/%s/%s/{pairprefix}_1%s" % (wildcards.datatype,
                                                                                                         wildcards.stage,
                                                                                                         config["fastq_extension"])),
                                                pairprefix=input_pairprefix_dict[wildcards.datatype]),
        reverse_fastqs=lambda wildcards: expand(output_dict["data"] / ("fastq/%s/%s/{pairprefix}_2%s" % (wildcards.datatype,
                                                                                                         wildcards.stage,
                                                                                                         config["fastq_extension"])),
                                                pairprefix=input_pairprefix_dict[wildcards.datatype]),
    output:
        db=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.yak" # , (?!^histo$)
    log:
        std=output_dict["log"] / "yak.{datatype}.{stage}.{kmer_length}.log",
        cluster_log=output_dict["cluster_log"] / "yak.{datatype}.{stage}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "yak.{datatype}.{stage}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "yak.{datatype}.{stage}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("yak_pe_reads"),
        cpus=get_threads(parameters["threads"]["yak_pe_reads"], "cpu"),
        time=parameters["time"]["yak_pe_reads"],
        mem=lambda wildcards, attempt: attempt * parameters["memory_mb"]["yak_pe_reads"],
        yak_kmer_counter=1
    threads:
        parameters["threads"]["yak_pe_reads"]
    shell:
         " yak count -o {output} -k {wildcards.kmer_length} -b 37 <(zcat {input.forward_fastqs}) <(zcat {input.reverse_fastqs}) > {log.std} 2>&1; "
