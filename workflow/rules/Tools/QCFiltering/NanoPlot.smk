localrules: gather_nanoplot_stats_per_stage, gather_nanoplot_stats_per_stage_track_data

rule nanoplot:
    input:
        fastq=config["out_dir"] / ("data/{longread_datatype}/{stage}/{fileprefix}%s" % config["fastq_ext"])
    output:
        yield_png=config["out_dir"] / "qc/nanoplot/{longread_datatype}/{stage}/{fileprefix}.Yield_By_Length.png",
        stats=config["out_dir"] / "qc/nanoplot/{longread_datatype}/{stage}/{fileprefix}.NanoStats.txt",
        pickle=config["out_dir"] / "qc/nanoplot/{longread_datatype}/{stage}/{fileprefix}.NanoPlot-data.pickle"
    log:
        std=config["out_dir"] / "log/nanoplot.{longread_datatype}.{stage}.{fileprefix}.log",
        cluster_log=config["out_dir"] / "log/nanoplot.{longread_datatype}.{stage}.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/nanoplot.{longread_datatype}.{stage}.{fileprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/nanoplot.{longread_datatype}.{stage}.{fileprefix}.benchmark.txt"
    conda:
        config["conda"]["nanopore"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["nanopore"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("nanoplot"),
        cpus=parameters["threads"]["nanoplot"],
        time=parameters["time"]["nanoplot"],
        mem=parameters["memory_mb"]["nanoplot"],
    threads:
        parameters["threads"]["nanoplot"]
    shell:
        " NanoPlot -f png svg -t {threads} --store --tsv_stats  -o `dirname {output.stats}` -p {wildcards.fileprefix}. "
        " --plots kde dot  --dpi 300 --fastq {input.fastq} > {log.std} 2>&1; "

use rule nanoplot as nanoplot_track_data with:
    input:
        fastq=config["out_dir"] / ("track_data/{longread_datatype}/{track_name}/{stage}/{fileprefix}%s" % config["fastq_ext"])
    output:
        yield_png=config["out_dir"] / "track_qc/nanoplot/{longread_datatype}/{track_name}/{stage}/{fileprefix}.Yield_By_Length.png",
        stats=config["out_dir"] / "track_qc/nanoplot/{longread_datatype}/{track_name}/{stage}/{fileprefix}.NanoStats.txt",
        pickle=config["out_dir"] / "track_qc/nanoplot/{longread_datatype}/{track_name}/{stage}/{fileprefix}.NanoPlot-data.pickle"
    log:
        std=config["out_dir"] / "log/nanoplot_track_data.{longread_datatype}.{track_name}.{stage}.{fileprefix}.log",
        cluster_log=config["out_dir"] / "log/nanoplot_track_data.{longread_datatype}.{track_name}.{stage}.{fileprefix}.cluster.log",
        cluster_err=config["out_dir"] / "log/nanoplot_track_data.{longread_datatype}.{track_name}.{stage}.{fileprefix}.cluster.err"
    benchmark:
        config["out_dir"] / "log/nanoplot_track_data.{longread_datatype}.{track_name}.{stage}.{fileprefix}.benchmark.txt"

rule gather_nanoplot_stats_per_stage:
    input:
        stats=lambda wildcards: expand(rules.nanoplot.output.stats,
                                       fileprefix=config["data"][wildcards.longread_datatype]["conv_file_prefix_list"], allow_missing=True)
    output:
        stage_stats=config["out_dir"] / "qc/nanoplot/{longread_datatype}/{stage}/{longread_datatype}.{stage}.NanoStats.tsv",
    params:
        labels=lambda wildcards: ",".join(config["data"][wildcards.longread_datatype]["file_prefix_list"])
    log:
        std=config["out_dir"] / "log/gather_nanoplot_stats_per_stage.{longread_datatype}.{stage}.log",
        cluster_log=config["out_dir"] / "log/gather_nanoplot_stats_per_stage.{longread_datatype}.{stage}.cluster.log",
        cluster_err=config["out_dir"] / "log/gather_nanoplot_stats_per_stage.{longread_datatype}.{stage}.cluster.err"
    benchmark:
        config["out_dir"] / "log/gather_nanoplot_stats_per_stage.{longread_datatype}.{stage}.benchmark.txt"
    conda:
        config["conda"]["nanopore"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["nanopore"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("nanoplot"),
        cpus=parameters["threads"]["nanoplot"],
        time=parameters["time"]["nanoplot"],
        mem=parameters["memory_mb"]["nanoplot"],
    threads:
        parameters["threads"]["nanoplot"]
    shell:
        " workflow/scripts/stats/gather_nanoplot_output.py -s {wildcards.stage} -d {wildcards.longread_datatype} "
        "     -l {params.labels} -g -o {output.stage_stats} {input.stats} > {log.std} 2>&1; "

use rule gather_nanoplot_stats_per_stage as gather_nanoplot_stats_per_stage_track_data with:
    input:
        stats=lambda wildcards: expand(rules.nanoplot_track_data.output.stats,
                                       fileprefix=config["track_data"][wildcards.longread_datatype][wildcards.track_name]["conv_file_prefix_list"], allow_missing=True)
    output:
        stage_stats=config["out_dir"] / "track_qc/nanoplot/{longread_datatype}/{track_name}/{stage}/{longread_datatype}.{track_name}.{stage}.NanoStats.tsv",
    params:
        labels=lambda wildcards: ",".join(config["track_data"][wildcards.longread_datatype][wildcards.track_name]["file_prefix_list"])
    log:
        std=config["out_dir"] / "log/gather_nanoplot_stats_per_stage_track_data.{longread_datatype}.{track_name}.{stage}.log",
        cluster_log=config["out_dir"] / "log/gather_nanoplot_stats_per_stage_track_data.{longread_datatype}.{track_name}.{stage}.cluster.log",
        cluster_err=config["out_dir"] / "log/gather_nanoplot_stats_per_stage_track_data.{longread_datatype}.{track_name}.{stage}.cluster.err"
    benchmark:
        config["out_dir"] / "log/gather_nanoplot_stats_per_stage_track_data.{longread_datatype}.{track_name}.{stage}.benchmark.txt"