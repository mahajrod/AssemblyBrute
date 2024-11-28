localrules: gather_nanoplot_stats_per_stage

rule nanoplot:
    input:
        fastq=output_dict["data"] / ("fastq/{datatype}/{stage}/{fileprefix}%s" % config["fastq_extension"])
    output:
        yield_png=output_dict["qc"] / "nanoplot/{datatype, [^/]+}/{stage, [^/]+}/{fileprefix, [^/]+}.Yield_By_Length.png",
        stats=output_dict["qc"] / "nanoplot/{datatype, [^/]+}/{stage, [^/]+}/{fileprefix, [^/]+}.NanoStats.txt"
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

rule gather_nanoplot_stats_per_stage:
    input:
        stats=lambda wildcards: expand(rules.nanoplot.output.stats,
                                       fileprefix=input_file_prefix_dict[wildcards.datatype], allow_missing=True)
    output:
        stage_stats=output_dict["qc"] / "nanoplot/{datatype, [^/]+}/{stage, [^/]+}/{datatype}.{stage}.NanoStats.tsv",
    params:
        labels=lambda wildcards: ",".join(input_file_prefix_dict[wildcards.datatype])
    log:
        std=output_dict["log"]/ "nanoplot.{datatype}.{stage}.log",
        #stats=log_dir_path / "{library_id}/fastqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"] / "nanoplot.{datatype}.{stage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "nanoplot.{datatype}.{stage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "nanoplot.{datatype}.{stage}.benchmark.txt"
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
        " workflow/scripts/stats/gather_nanoplot_output.py -s {wildcards.stage} -d {wildcards.datatype} "
        " -l {params.labels} -g -o {output.stage_stats} {input.stats} > {log.std} 2>&1; "

"""
rule gather_datatype_nanoplot_stats:
    input:
        stats=lambda wildcards: expand(rules.gather_nanoplot_stats_per_stage.output.stage_stats,
                                       datatype_list=long_read_data_type_set,
                                       stage_list=
                                      )
    output:
        stage_stats=output_dict["qc"] / "nanoplot/{datatype, [^/]+}/{stage, [^/]+}/NanoStats.tsv",
    params:
        labels=lambda wildcards: ",".join(input_file_prefix_dict[wildcards.datatype])
    log:
        std=output_dict["log"]/ "nanoplot.{datatype}.{stage}.log",
        #stats=log_dir_path / "{library_id}/fastqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"] / "nanoplot.{datatype}.{stage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "nanoplot.{datatype}.{stage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "nanoplot.{datatype}.{stage}.benchmark.txt"
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
        " workflow/scripts/stats/gather_nanoplot_output.py -s {wildcards.stage} -d {wildcards.datatype} "
        " -l {params.labels} -g -o {output.stage_stats} {input.stats} > {log.std} 2>&1; "

"""