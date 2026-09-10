
rule gfa2og:
    input:
        gfa="{gfa_dir}/{gfa_prefix}.gfa",
        log_dir=ancient("{gfa_dir}/log/")
    output:
        og="{gfa_dir}/{gfa_prefix}.og"
    log:
        std="{gfa_dir}/log/gfa2og.{gfa_prefix}.log",
        cluster_log="{gfa_dir}/log/gfa2og.{gfa_prefix}.cluster.log",
        cluster_err="{gfa_dir}/log/gfa2og.{gfa_prefix}.cluster.err"
    benchmark:
        "{gfa_dir}/log/gfa2og.{gfa_prefix}.benchmark.txt"
    conda:
        config["conda"]["hifiasm"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["hifiasm"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("gfa2og"),
        cpus=parameters["threads"]["gfa2og"],
        time=parameters["time"]["gfa2og"],
        mem=parameters["memory_mb"]["gfa2og"],
    threads:
        parameters["threads"]["gfa2og"]
    shell:
         " odgi build -t {threads} -g {input.gfa} -o {output.og} >{log.std} 2>&1; "

rule og_sort_1d:
    input:
        og="{gfa_dir}/{gfa_prefix}.og",
        log_dir=ancient("{gfa_dir}/log/")
    output:
        sorted_1d_og="{gfa_dir}/{gfa_prefix}.sorted_1d.og"
    log:
        std="{gfa_dir}/log/og_sort_1d.{gfa_prefix}.log",
        cluster_log="{gfa_dir}/log/og_sort_1d.{gfa_prefix}.cluster.log",
        cluster_err="{gfa_dir}/log/og_sort_1d.{gfa_prefix}.cluster.err"
    benchmark:
        "{gfa_dir}/log/og_sort_1d{gfa_prefix}.benchmark.txt"
    conda:
        config["conda"]["hifiasm"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["hifiasm"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("og_sort_1d"),
        cpus=parameters["threads"]["og_sort_1d"],
        time=parameters["time"]["og_sort_1d"],
        mem=parameters["memory_mb"]["og_sort_1d"],
    threads:
        parameters["threads"]["og_sort_1d"]
    shell:
         " odgi sort -i {input.og} --threads {threads} -P -Y -o {output.sorted_1d_og} > {log.std} 2>&1; "

rule og_viz_1d:
    input:
        sorted_1d_og="{gfa_dir}/{gfa_prefix}.sorted_1d.og",
        log_dir=ancient("{gfa_dir}/log/")
    output:
        sorted_1d_png="{gfa_dir}/{gfa_prefix}.sorted_1d.png",
    log:
        std="{gfa_dir}/log/og_viz_1d.{gfa_prefix}.log",
        cluster_log="{gfa_dir}/log/og_viz_1d.{gfa_prefix}.cluster.log",
        cluster_err="{gfa_dir}/log/og_viz_1d.{gfa_prefix}.cluster.err"
    benchmark:
        "{gfa_dir}/log/og_viz_1d{gfa_prefix}.benchmark.txt"
    conda:
        config["conda"]["hifiasm"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["hifiasm"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("og_viz_1d"),
        cpus=parameters["threads"]["og_viz_1d"],
        time=parameters["time"]["og_viz_1d"],
        mem=parameters["memory_mb"]["og_viz_1d"],
    threads:
        parameters["threads"]["og_viz_1d"]
    shell:
         " odgi viz --threads {threads} -P -i {input.sorted_1d_og} -o {output.sorted_1d_png} > {log.std} 2>&1; "
