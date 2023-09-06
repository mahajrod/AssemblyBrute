
rule fastqc:
    input:
        #fastq_dir=rules.create_fastq_links.output,
        fastq=output_dict["data"] / ("fastq/{datatype}/{stage}/{fileprefix}%s" % config["fastq_extension"])
    output:
        zip=output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip" ,
        #stats=merged_raw_fastqc_dir_path / "{library_id}/{library_id}.raw.fastqc.stats"
    params:
        kmer=lambda wildcards: parse_option("kmer_length", parameters["tool_options"]["fastqc"][wildcards.datatype], " -k "),
        #out_dir=lambda wildcards: output_dict["qc"] / "fastqc/{0}/{1}/".format(wildcards.datatype,
        #                                                                       wildcards.stage),
        nogroup=lambda wildcards: "" if wildcards.datatype in config["long_read_data"] else "--nogroup" # turns off base grouping for short reads

    log:
        std=output_dict["log"]/ "fastqc.{datatype}.{stage}.{fileprefix}.log",
        #stats=log_dir_path / "{library_id}/fastqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"] / "fastqc.{datatype}.{stage}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "fastqc.{datatype}.{stage}.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "fastqc.{datatype}.{stage}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("fastqc"),
        cpus=parameters["threads"]["fastqc"],
        time=parameters["time"]["fastqc"],
        mem=parameters["memory_mb"]["fastqc"],
    threads:
        parameters["threads"]["fastqc"]
    shell:
        " OUTDIR=`dirname {output.zip}`;"
        " fastqc {params.nogroup} --memory {resources.mem} {params.kmer} -t {threads} "
        " -o ${{OUTDIR}} {input} 1>{log.std} 2>&1; "
