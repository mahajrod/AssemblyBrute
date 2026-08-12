
rule nanoqc:
    input:
        fastq=config["out_dir"] / ("data/{longread_datatype}/{stage}/{fileprefix}%s" % config["fastq_ext"])
    output:
        dir=directory(config["out_dir"] / "qc/nanoqc/{longread_datatype}/{stage}/{fileprefix}")
    log:
        std=config["out_dir"] / "log/nanoqc.{longread_datatype}.{stage}.{fileprefix}.log",
        cluster_log=config["out_dir"] / "log/nanoqc.{longread_datatype}.{stage}.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/nanoqc.{longread_datatype}.{stage}.{fileprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/nanoqc.{longread_datatype}.{stage}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["nanopore"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["nanopore"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("nanoqc"),
        cpus=parameters["threads"]["nanoqc"],
        time=parameters["time"]["nanoqc"],
        mem=parameters["memory_mb"]["nanoqc"],
    threads:
        parameters["threads"]["nanoqc"]
    shell:
        "mkdir -p {output.dir};  nanoQC -o {output.dir} {input} 1>{log.std} 2>&1; "

use rule nanoqc as nanoqc_track_data with:
    input:
        fastq=config["out_dir"] / ("ext_track_data/{longread_datatype}/{track_name}/{stage}/{fileprefix}%s" % config["fastq_ext"])
    output:
        dir=directory(config["out_dir"] / "ext_track_qc/nanoqc/{longread_datatype}/{track_name}/{stage}/{fileprefix}")
    log:
        std=config["out_dir"] / "log/nanoqc_track_data.{longread_datatype}.{track_name}.{stage}.{fileprefix}.log",
        cluster_log=config["out_dir"] / "log/nanoqc_track_data.{longread_datatype}.{track_name}.{stage}.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/nanoqc_track_data.{longread_datatype}.{track_name}.{stage}.{fileprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/nanoqc_track_data.{longread_datatype}.{track_name}.{stage}{fileprefix}.benchmark.txt"