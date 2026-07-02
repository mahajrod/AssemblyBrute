
def get_files_for_jellyfish(wildcards):
    file_list = []
    for datatype in wildcards.datatype.split("_"):
        file_list += expand(config["out_dir"] / ("data/%s/%s/{fileprefix}%s" % (datatype,
                                                                                          wildcards.stage,
                                                                                          config["data"][wildcards.datatype]["conv_ext"])),
                            fileprefix=config["data"][wildcards.datatype]["conv_file_prefix_list"],
                            allow_missing=True)

    return file_list

rule jellyfish:
    input:get_files_for_jellyfish

    output:
        jf=config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.jellyfish.jf",
        histo=config["out_dir"] / "kmer/{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.jellyfish.histo"
    params:
        hash_size=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["hash_size"],
        min_coverage=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["min_coverage"],
        max_coverage=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["max_coverage"],
        increment=lambda wildcards: parameters["tool_options"]["jellyfish"][wildcards.datatype]["increment"]
    log:
        count_log=config["out_dir"] / "log/jellyfish.{datatype}.{stage}.{kmer_length}.count.log",
        histo_log=config["out_dir"] / "log/jellyfish.{datatype}.{stage}.{kmer_length}.histo.log",
        cluster_log=config["out_dir"] / "log/jellyfish.{datatype}.{stage}.{kmer_length}.cluster.log",
        cluster_err=config["out_dir"] / "log/jellyfish.{datatype}.{stage}.{kmer_length}.cluster.err"
    benchmark:
        config["out_dir"] / "jellyfish.{datatype}.{stage}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("jellyfish"),
        cpus=get_threads(parameters["threads"]["jellyfish"], "cpu"),
        time=parameters["time"]["jellyfish"],
        mem=parameters["memory_mb"]["jellyfish"],
        kmer_counter=1
    threads:
        parameters["threads"]["jellyfish"]
    shell:
         " jellyfish count -C -m {wildcards.kmer_length} -s {params.hash_size} -t {threads} -o {output.jf}  "
         "     <(zcat {input}) 1>{log.count_log} 2>&1; "
         " jellyfish histo -o {output.histo} -t {threads} -l {params.min_coverage} -h {params.max_coverage} "
         "     -i {params.increment} {output.jf} 1>{log.histo_log} 2>&1; "

