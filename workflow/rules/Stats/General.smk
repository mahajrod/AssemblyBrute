localrules: gather_stats_per_stage_parameter, gather_stage_stats

rule gather_stats_per_stage_parameter:
    input:
        summary=lambda wildcards: expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.busco5.{busco_lineage}.summary",
                                         busco_lineage=config["busco_lineage_list"],
                                         haplotype=stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["haplotype_list"],
                                         allow_missing=True) if not config["skip_busco"] else [],
        quast_dirs=lambda wildcards: expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/quast/{genome_prefix}.{assembly_stage}.{haplotype}",
                                            haplotype=stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["haplotype_list"],
                                            allow_missing=True),
        qv_file=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/merqury/{genome_prefix}.{assembly_stage}.qv",
        completeness_stats_file=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/merqury/{genome_prefix}.{assembly_stage}.completeness.stats",
    params:
        #busco_list=lambda wildcards: (" -b " + ",".join(expand(out_dir_path / ("%s/%s/assembly_qc/busco5/%s.%s.{haplotype}.busco5.{busco_lineage}.summary" % (wildcards.assembly_stage,
        #                                                                                                                                                      wildcards.parameters,
        #                                                                                                                                                      wildcards.genome_prefix,
        #                                                                                                                                                      wildcards.assembly_stage)),
        #               busco_lineage=config["busco_lineage_list"],
        #               haplotype=haplotype_list,
        #               allow_missing=True) )) if not config["skip_busco"] else "",
        haplotype_list=lambda wildcards: ",".join(stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["haplotype_list"]),
        busco_lineage_list=(" -b " + ",".join(config["busco_lineage_list"])) if not config["skip_busco"] else ""
    output:
        stats=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/{genome_prefix}.{assembly_stage}.parameter_stats"
    log:
        std=output_dict["log"]/ "gather_stats_per_stage_parameter.{genome_prefix}.{assembly_stage}.{parameters}.log",
        #stats=log_dir_path / "{library_id}/multiqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"]/ "gather_stats_per_stage_parameter{genome_prefix}.{assembly_stage}.{parameters}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "gather_stats_per_stage_parameter.{genome_prefix}.{assembly_stage}.{parameters}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "gather_stats_per_stage_parameter.{genome_prefix}.{assembly_stage}.{parameters}.benchmark.txt"
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
        " ./workflow/scripts/gather_qc_stats.py -q results/{wildcards.assembly_stage}/{wildcards.parameters}/assembly_qc/ "
        " -p {wildcards.parameters} -e {wildcards.genome_prefix}.{wildcards.assembly_stage} -s {wildcards.assembly_stage} "
        " -a {params.haplotype_list} {params.busco_lineage_list} -o {output.stats} > {log.std} 2>&1; "


rule gather_stage_stats:
    input:
        stats=lambda wildcards: expand(out_dir_path / ("%s/{parameters}/assembly_qc/%s.%s.parameter_stats" % (wildcards.assembly_stage,
                                                                                                    wildcards.genome_prefix,
                                                                                                    wildcards.assembly_stage)),
                                       parameters=stage_dict[wildcards.assembly_stage]["parameters"].keys(),
                                       allow_missing=True)

    output:
        stats=out_dir_path / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats"
    log:
        std=output_dict["log"]/ "gather_stage_stats.{genome_prefix}.{assembly_stage}.log",
        #stats=log_dir_path / "{library_id}/multiqc_merged_raw.stats.log",
        cluster_log=output_dict["cluster_log"]/ "gather_stage_stats.{genome_prefix}.{assembly_stage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "gather_stage_stats.{genome_prefix}.{assembly_stage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "gather_stage_stats.{genome_prefix}.{assembly_stage}.benchmark.txt"
    #conda:
    #    config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("gather_stage_stats"),
        cpus=parameters["threads"]["gather_stage_stats"],
        time=parameters["time"]["gather_stage_stats"],
        mem=parameters["memory_mb"]["gather_stage_stats"],
    threads:
        parameters["threads"]["gather_stage_stats"]
    run:
        df_list = [pd.read_csv(filename, sep="\t", header=0,) for filename in input.stats]
        #print(df_list)
        merged_df = pd.concat(df_list)
        columns = list(merged_df.columns)
        merged_df = merged_df[columns[1:3] + [columns[0]] + columns[3:]].sort_values(by=["stage", "parameters", "haplotype"])

        merged_df.to_csv(output.stats, sep="\t", header=True, index=False)
