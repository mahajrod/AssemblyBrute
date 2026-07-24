localrules: gather_stats_per_stage_parameter, gather_stage_stats

rule gather_stats_per_stage_parameter:
    input:
        summary=lambda wildcards: expand(config["out_dir"]  / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}.busco5.summary",
                                         busco_lineage=config["busco_lineage_list"],
                                         haplotype=stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"],
                                         allow_missing=True) if (not config["skip_busco"]) and (config["assembly_qc_level"][wildcards.assembly_stage] >= 5) else [],
        quast_dirs=lambda wildcards: expand(config["out_dir"]  / "{assembly_stage}/{parameters}/assembly_qc/quast/{genome_prefix}.{assembly_stage}.{haplotype}",
                                            haplotype=stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"],
                                            allow_missing=True),
        qv_files=lambda wildcards: expand(config["out_dir"]  / ("%s/%s/assembly_qc/merqury/{merqury_datatype}/%s.%s.{merqury_datatype}.qv" % (wildcards.assembly_stage,
                                                                                                                                                    wildcards.parameters,
                                                                                                                                                    wildcards.genome_prefix,
                                                                                                                                                    wildcards.assembly_stage)),
                                          merqury_datatype=parameters["tool_options"]["assembly_qc"]["merqury"]["datatype_list"],
                                          allow_missing=True) if config["assembly_qc_level"][wildcards.assembly_stage] >= 2 else [],
        completeness_stats_files=lambda wildcards: expand(config["out_dir"]  / ("%s/%s/assembly_qc/merqury/{merqury_datatype}/%s.%s.{merqury_datatype}.completeness.stats" % (wildcards.assembly_stage,
                                                                                                                                                    wildcards.parameters,
                                                                                                                                                    wildcards.genome_prefix,
                                                                                                                                                    wildcards.assembly_stage)),
                                                          merqury_datatype=parameters["tool_options"]["assembly_qc"]["merqury"]["datatype_list"],
                                                          allow_missing=True) if config["assembly_qc_level"][wildcards.assembly_stage] >= 2 else [],
    params:
        merqury_datatypes=lambda  wildcards: ("-m " + ",".join(parameters["tool_options"]["assembly_qc"]["merqury"]["datatype_list"])) if config["assembly_qc_level"][wildcards.assembly_stage] >= 1 else "",
        assembly_prefix_list=lambda wildcards: ",".join(map(lambda haplotype: f"{wildcards.genome_prefix}.{wildcards.assembly_stage}.{haplotype}",
                                                            stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"])),
        busco_lineage_list=lambda  wildcards: (" -b " + ",".join(config["busco_lineage_list"])) if (not config["skip_busco"]) and (config["assembly_qc_level"][wildcards.assembly_stage] >= 4) else ""
    output:
        stats=config["out_dir"]  / "{assembly_stage}/{parameters}/assembly_qc/{genome_prefix}.{assembly_stage}.parameter_stats"
    log:
        std=config["out_dir"] / "log/gather_stats_per_stage_parameter.{genome_prefix}.{assembly_stage}.{parameters}.log",
        cluster_log=config["out_dir"] / "log/gather_stats_per_stage_parameter{genome_prefix}.{assembly_stage}.{parameters}.cluster.log",
        cluster_err=config["out_dir"] / "log/gather_stats_per_stage_parameter.{genome_prefix}.{assembly_stage}.{parameters}.cluster.err"
    benchmark:
        config["out_dir"] / "log/gather_stats_per_stage_parameter.{genome_prefix}.{assembly_stage}.{parameters}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("gather_stage_stats"),
        cpus=parameters["threads"]["gather_stage_stats"],
        time=parameters["time"]["gather_stage_stats"],
        mem=parameters["memory_mb"]["gather_stage_stats"],
    threads:
        parameters["threads"]["gather_stage_stats"]
    shell:
        " ./workflow/scripts/gather_qc_stats.py -q results/{wildcards.assembly_stage}/{wildcards.parameters}/assembly_qc/ "
        "     -p {wildcards.parameters} -e {wildcards.genome_prefix}.{wildcards.assembly_stage} -s {wildcards.assembly_stage} "
        "     -a {params.assembly_prefix_list} {params.busco_lineage_list} -o {output.stats} {params.merqury_datatypes} > {log.std} 2>&1; "


rule gather_stage_stats:
    priority: 100000
    input:
        stats=lambda wildcards: expand(config["out_dir"]  / ("%s/{parameters}/assembly_qc/%s.%s.parameter_stats" % (wildcards.assembly_stage,
                                                                                                    wildcards.genome_prefix,
                                                                                                    wildcards.assembly_stage)),
                                       parameters=stage_dict[wildcards.assembly_stage].parameters.keys(),
                                       allow_missing=True)

    output:
        stats=config["out_dir"] / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats"
    log:
        std=config["out_dir"] / "log/gather_stage_stats.{genome_prefix}.{assembly_stage}.log",
        cluster_log=config["out_dir"] / "log/gather_stage_stats.{genome_prefix}.{assembly_stage}.cluster.log",
        cluster_err=config["out_dir"] / "log/gather_stage_stats.{genome_prefix}.{assembly_stage}.cluster.err"
    benchmark:
        config["out_dir"] / "log/gather_stage_stats.{genome_prefix}.{assembly_stage}.benchmark.txt"
    #conda:
    #    config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("gather_stage_stats"),
        cpus=parameters["threads"]["gather_stage_stats"],
        time=parameters["time"]["gather_stage_stats"],
        mem=parameters["memory_mb"]["gather_stage_stats"],
    threads:
        parameters["threads"]["gather_stage_stats"]
    run:
        df_list = [pd.read_csv(filename, sep="\t", header=0,) for filename in input.stats]
        merged_df = pd.concat(df_list)
        columns = list(merged_df.columns)
        merged_df = merged_df[columns[1:3] + [columns[0]] + columns[3:]].sort_values(by=["stage", "parameters", "assembly_prefix"])

        merged_df.to_csv(output.stats, sep="\t", header=True, index=False)
