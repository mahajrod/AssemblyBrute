
rule nanoqc:
    input:
        fastq=output_dict["data"] / ("fastq/{datatype}/{stage}/{fileprefix}%s" % config["fastq_extension"])
    output:
        dir=directory(output_dict["qc"] / "nanoqc/{datatype}/{stage}/{fileprefix}")
    log:
        std=output_dict["log"]/ "nanoqc.{datatype}.{stage}.{fileprefix}.log",
        #stats=log_dir_path / "{library_id}/fastqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"] / "nanoqc.{datatype}.{stage}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "nanoqc.{datatype}.{stage}.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "nanoqc.{datatype}.{stage}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["nanopore"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["nanopore"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("nanoqc"),
        cpus=parameters["threads"]["nanoqc"],
        time=parameters["time"]["nanoqc"],
        mem=parameters["memory_mb"]["nanoqc"],
    threads:
        parameters["threads"]["nanoqc"]
    shell:
        "mkdir -p {output.dir};  nanoQC -o {output.dir} {input} 1>{log.std} 2>&1; "