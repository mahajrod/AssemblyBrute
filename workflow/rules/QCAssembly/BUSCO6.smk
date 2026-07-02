
rule busco6_download: #
    priority: 500
    output:
        lineage_dir=directory(out_dir_path / "download/busco5/lineages/{busco6_lineage}"),
    params:
        busco_download_dir=out_dir_path / "download/busco5/"
    log:
        std=output_dict["log"] / "busco5_download.{busco6_lineage}.log",
        cluster_log=output_dict["cluster_log"] / "busco5_download.{busco6_lineage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "busco5_download.{busco6_lineage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "busco6_download.{busco6_lineage}.benchmark.txt"
    conda:
        config["conda"]["busco6"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco6"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("busco_download"),
        cpus=parameters["threads"]["busco_download"],
        time=parameters["time"]["busco_download"],
        mem=parameters["memory_mb"]["busco_download"],
    threads:
        parameters["threads"]["busco_download"]
    shell:
         " busco --download_path {params.busco_download_dir} --download {wildcards.busco6_lineage} > {log.std} 2>&1; "
