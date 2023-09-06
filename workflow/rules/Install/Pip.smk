localrules: install_from_pip

rule install_from_pip:
    priority: 100000
    input:
        requirements="workflow/envs/pip.{conda_env}.requirements"
    output:
        requirements=temp("results/config/pip.{conda_env}.requirements") # temp to force rule execution every... hm, nearly every... run of the pipeline
    log:
        pip=output_dict["log"]  / "install_from_pip.{conda_env}.pip.log",
        cp=output_dict["log"]  / "install_from_pip.{conda_env}.cp.log",
        cluster_log=output_dict["cluster_log"] / "install_from_pip.{conda_env}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "install_from_pip.{conda_env}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "install_from_pip.{conda_env}.benchmark.txt"
    conda:
        lambda wildcards: config["conda"][wildcards.conda_env]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"][wildcards.conda_env]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("install_from_pip"),
        cpus=parameters["threads"]["install_from_pip"] ,
        time=parameters["time"]["install_from_pip"],
        mem=parameters["memory_mb"]["install_from_pip"]
    threads: parameters["threads"]["install_from_pip"]

    shell:
        " pip install -r {input.requirements} > {log.pip} 2>&1 ;"
        " cp -f {input.requirements} {output.requirements} > {log.cp} 2>&1;"