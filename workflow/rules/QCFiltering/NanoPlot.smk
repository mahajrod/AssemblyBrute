
rule nanoplot:
    input:
        fastq=output_dict["data"] / ("fastq/{datatype}/{stage}/{fileprefix}%s" % config["fastq_extension"])
    output:
        yield_png=output_dict["qc"] / "nanoplot/{datatype}/{stage}/{fileprefix}.Yield_By_Length.png"

    log:
        std=output_dict["log"]/ "nanoplot.{datatype}.{stage}.{fileprefix}.log",
        #stats=log_dir_path / "{library_id}/fastqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"] / "nanoplot.{datatype}.{stage}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "nanoplot.{datatype}.{stage}.{fileprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "nanoplot.{datatype}.{stage}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["nanopore"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["nanopore"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("nanoplot"),
        cpus=parameters["threads"]["nanoplot"],
        time=parameters["time"]["nanoplot"],
        mem=parameters["memory_mb"]["nanoplot"],
    threads:
        parameters["threads"]["nanoplot"]
    shell:
        " NanoPlot -f png svg -t {threads} --store --tsv_stats  -o `dirname {output}` -p {wildcards.fileprefix}. "
        " --plots kde dot  --dpi 300 --fastq {input.fastq} > {log.std} 2>&1; "
