#DEPRECATED. Will be removed soon
rule kat_gcp:
    input:
        jf=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.jellyfish.jf"
    output:
        png=output_dict["kmer"] / "{datatype}/{stage}/kat/{datatype}.{stage}.{kmer_length}.jellyfish.kat.gcp.mx.png",
    params:
        #kmer_length=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["kmer_length"],
        hash_size=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["hash_size"],
        min_coverage=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["min_coverage"],
        max_coverage=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["max_coverage"],
        increment=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["increment"]
    log:
        std=output_dict["log"] / "jellyfish.{datatype}.{stage}.{kmer_length}.log",
        cluster_log=output_dict["cluster_log"] / "jellyfish.{datatype}.{stage}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "jellyfish.{datatype}.{stage}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "jellyfish.{datatype}.{stage}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["kat"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kat"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        cpus=parameters["threads"]["kat_gcp"],
        time=parameters["time"]["kat_gcp"],
        mem=parameters["memory_mb"]["kat_gcp"],
    threads:
        parameters["threads"]["kat_gcp"]
    shell:
         " PNG={output.png}; "
         " kat gcp -t {threads} -m {wildcards.kmer_length} -o ${{PNG%.mx.png}} {input.jf}  1>{log.std} 2>&1;"

