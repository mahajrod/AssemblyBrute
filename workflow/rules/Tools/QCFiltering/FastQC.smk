
rule fastqc:
    input:
        fastq=config["out_dir"] / ("data/{fastqc_datatype}/{stage}/{fileprefix}%s" % config["fastq_ext"])
    output:
        zip=config["out_dir"] / "qc/fastqc/{fastqc_datatype}/{stage}/{fileprefix}_fastqc.zip"
    params:
        kmer=lambda wildcards: parse_option("kmer_length", parameters["tool_options"]["fastqc"][wildcards.fastqc_datatype], " -k "),
        nogroup=lambda wildcards: "" if wildcards.fastqc_datatype in config["data_feature_dict"]["long_read"] else "--nogroup" # turns off base grouping for short reads

    log:
        std=config["out_dir"] / "log/fastqc.{fastqc_datatype}.{stage}.{fileprefix}.log",
        cluster_log=config["out_dir"] / "log/fastqc.{fastqc_datatype}.{stage}.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/fastqc.{fastqc_datatype}.{stage}.{fileprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/fastqc.{fastqc_datatype}.{stage}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("fastqc"),
        cpus=parameters["threads"]["fastqc"],
        time=parameters["time"]["fastqc"],
        mem=parameters["memory_mb"]["fastqc"],
    threads:
        parameters["threads"]["fastqc"]
    shell:
        " OUTDIR=`dirname {output.zip}`;"
        " fastqc {params.nogroup} --memory {resources.mem} {params.kmer} -t {threads} "
        "     -o ${{OUTDIR}} {input.fastq} 1>{log.std} 2>&1; "

use rule fastqc as fastqc_track_data with:
    input:
        fastq=config["out_dir"] / ("ext_track_data/{fastqc_datatype}/{track_name}/{stage}/{fileprefix}%s" % config["fastq_ext"])
    output:
        zip=config["out_dir"] / "ext_track_qc/fastqc/{fastqc_datatype}/{track_name}/{stage}/{fileprefix}_fastqc.zip"
    log:
        std=config["out_dir"] / "log/fastqc_track_data.{fastqc_datatype}.{track_name}.{stage}.{fileprefix}.log",
        cluster_log=config["out_dir"] / "log/fastqc_track_data.{fastqc_datatype}.{track_name}.{stage}.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/fastqc_track_data.{fastqc_datatype}.{track_name}.{stage}.{fileprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/fastqc_track_data.{fastqc_datatype}.{track_name}.{stage}.{fileprefix}.benchmark.txt"


