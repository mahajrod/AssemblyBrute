rule jellyfish:
    input:
        lambda wildcards:
            expand(output_dict["data"] / ("fastq/%s/%s/{fileprefix}%s" % (wildcards.datatype,
                                                                          wildcards.stage,
                                                                          config["fastq_extension"])),
                   fileprefix=input_file_prefix_dict[wildcards.datatype],
                   allow_missing=True
                   )
    output:
        jf=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.jellyfish.jf",
        histo=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.jellyfish.histo"
    params:
        #kmer_length=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["kmer_length"],
        hash_size=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["hash_size"],
        min_coverage=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["min_coverage"],
        max_coverage=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["max_coverage"],
        increment=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["increment"]
    log:
        count_log=output_dict["log"] / "jellyfish.{datatype}.{stage}.{kmer_length}.count.log",
        histo_log=output_dict["log"] / "jellyfish.{datatype}.{stage}.{kmer_length}.histo.log",
        cluster_log=output_dict["cluster_log"] / "jellyfish.{datatype}.{stage}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "jellyfish.{datatype}.{stage}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "jellyfish.{datatype}.{stage}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["jellyfish"],
        time=parameters["time"]["jellyfish"],
        mem=parameters["memory_mb"]["jellyfish"],
        kmer_counter=1
    threads:
        parameters["threads"]["jellyfish"]
    shell:
         " jellyfish count -C -m {wildcards.kmer_length} -s {params.hash_size} -t {threads} -o {output.jf}  "
         " <(zcat {input}) 1>{log.count_log} 2>&1; "
         " jellyfish histo -o {output.histo} -t {threads} -l {params.min_coverage} -h {params.max_coverage} "
         " -i {params.increment} {output.jf} 1>{log.histo_log} 2>&1;"

