localrules: create_local_log_dir, create_local_benchmark_dir

rule create_local_log_dir:
    input:
        work_dir="{directory}"
    output:
        log_dir=directory("{directory}/log/")
    log:
        std="{directory}/log.log",
        cluster_log="{directory}/log.cluster.log",
        cluster_err="{directory}/log.cluster.err"
    benchmark:
        "{directory}/log.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("gather_stage_stats"),
        cpus=parameters["threads"]["gather_stage_stats"],
        time=parameters["time"]["gather_stage_stats"],
        mem=parameters["memory_mb"]["gather_stage_stats"],
    threads:
        parameters["threads"]["gather_stage_stats"]
    shell:
        " mkdir -p {output.log_dir} > {log.std} 2>&1; "

rule create_local_benchmark_dir:
    input:
        work_dir="{directory}",
        log_dir="{directory}/log/"
    output:
        benchmark_dir=directory("{directory}/benchmark/")
    log:
        std="{directory}/log/create_local_benchmark_dir.log",
        cluster_log="{directory}/log/create_local_benchmark_dir.log",
        cluster_err="{directory}/log/create_local_benchmark_dir.err"
    benchmark:
        "{directory}/log/create_local_benchmark_dir.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("gather_stage_stats"),
        cpus=parameters["threads"]["gather_stage_stats"],
        time=parameters["time"]["gather_stage_stats"],
        mem=parameters["memory_mb"]["gather_stage_stats"],
    threads:
        parameters["threads"]["gather_stage_stats"]
    shell:
        " mkdir -p {output.benchmark_dir} > {log.std} 2>&1; "
