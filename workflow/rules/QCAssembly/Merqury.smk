
#stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"]

rule merqury: # TODO: add handling for cases of haploid and polyploid genomes
    input:
        meryl_db_dir=output_dict["kmer"] / "{0}/{1}/{0}.{1}.{2}.meryl".format(config["final_kmer_datatype"],
                                                                              "filtered" if config["final_kmer_datatype"] in config["filtered_data"] else "raw",
                                                                              config["final_kmer_length"],) ,
        primary_assembly=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(wildcards.assembly_stage,
                                                                                             wildcards.parameters,
                                                                                             wildcards.genome_prefix,
                                                                                             stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["haplotype_list"][0]),
        alternative_assembly=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(wildcards.assembly_stage,
                                                                                             wildcards.parameters,
                                                                                             wildcards.genome_prefix,
                                                                                             "hap2") if stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["haplotype_list"][0] != "hap0" else [],
    output:
        qv_file=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/merqury/{genome_prefix}.{assembly_stage}.qv",
        completeness_stats_file=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/merqury/{genome_prefix}.{assembly_stage}.completeness.stats",
    params:
        #dir=lambda wildcards: out_dir_path / "{0}/{1}/assembly_qc/merqury/".format(wildcards.assembly_stage,
        #                                                                           wildcards.parameters),
        out_prefix=lambda wildcards: "{0}.{1}".format(wildcards.genome_prefix,
                                                      wildcards.assembly_stage)
    log:
        std=output_dict["log"].resolve() / "merqury.{assembly_stage}.{parameters}.{genome_prefix}.log",
        mkdir_log=(output_dict["log"]).resolve() / "merqury.{assembly_stage}.{parameters}.{genome_prefix}.mkdir.log",
        cd_log=(output_dict["log"]).resolve() / "merqury.{assembly_stage}.{parameters}.{genome_prefix}.cd.log",
        cluster_log=(output_dict["cluster_log"]).resolve() / "merqury.{assembly_stage}.{parameters}.{genome_prefix}.cluster.log",
        cluster_err=(output_dict["cluster_error"]).resolve() / "merqury.{assembly_stage}.{parameters}.{genome_prefix}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "merqury.{assembly_stage}.{parameters}.{genome_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("merqury"),
        cpus=parameters["threads"]["merqury"],
        time=parameters["time"]["merqury"],
        mem=parameters["memory_mb"]["merqury"],
    threads:
        parameters["threads"]["merqury"]
    shell:
         " OUT_DIR=`dirname {output.qv_file}`; "
         " MERYL_DB=`realpath -s {input.meryl_db_dir}`;"
         " PRIMARY_ASSEMBLY=`realpath -s {input.primary_assembly}`;"
         " if [ -z '{input.alternative_assembly}' ]; "
         " then "
         "      ALTERNATIVE_ASSEMBLY=''; "
         " else "
         "      ALTERNATIVE_ASSEMBLY=`realpath -s {input.alternative_assembly}`; "
         " fi;"
         " cd ${{OUT_DIR}}; "
         " OMP_NUM_THREADS={threads} merqury.sh ${{MERYL_DB}} "
         " ${{PRIMARY_ASSEMBLY}} ${{ALTERNATIVE_ASSEMBLY}} {params.out_prefix}  1>{log.std} 2>&1;"


