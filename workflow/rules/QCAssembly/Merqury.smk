
#stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"]

rule merqury: # TODO: add handling for cases of haploid and polyploid genomes
    input:
        meryl_db_dir=output_dict["kmer"] / "{0}/filtered/{0}.filtered.{1}.meryl".format(config["final_kmer_datatype"],
                                                                                        config["final_kmer_length"],),
        primary_assembly=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(wildcards.assembly_stage,
                                                                                             wildcards.parameters,
                                                                                             wildcards.genome_prefix,
                                                                                             "hap1" if stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["option_set"]["assembly_ploidy"] > 1 else "hap0"),
        alternative_assembly=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(wildcards.assembly_stage,
                                                                                             wildcards.parameters,
                                                                                             wildcards.genome_prefix,
                                                                                             "hap2") if stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["option_set"]["assembly_ploidy"] > 1 else [],
    output:
        qv_file=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/merqury/{genome_prefix}.{assembly_stage}.qv",
        completeness_stats_file=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/merqury/{genome_prefix}.{assembly_stage}.completeness.stats",
    params:
        #dir=lambda wildcards: out_dir_path / "{0}/{1}/assembly_qc/merqury/".format(wildcards.assembly_stage,
        #                                                                           wildcards.parameters),
        primary_assembly=lambda wildcards: f"`realpath -s {wildcards.primary_assembly}` ",
        #alternative_assembly=,
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
        cpus=parameters["threads"]["merqury"],
        time=parameters["time"]["merqury"],
        mem=parameters["memory_mb"]["merqury"],
    threads:
        parameters["threads"]["merqury"]
    shell:
         " OUT_DIR=`dirname {output.qv_file}`; "
         " MERYL_DB=`realpath -s {input.meryl_db_dir}`;"
         " if [ -z '{input.alternative_assembly}' ]; "
         " then "
         "      ALTERNATIVE_ASSEMBLY=''; "
         " else "
         "      ALTERNATIVE_ASSEMBLY=`realpath -s {input.alternative_assembly}`; "
         " fi"
         " PRIMARY_ASSEMBLY=`realpath -s {input.primary_assembly}`;"
         " cd ${{OUT_DIR}}; "
         " OMP_NUM_THREADS={threads} merqury.sh ${{MERYL_DB}} "
         " ${{PRIMARY_ASSEMBLY}} ${{ALTERNATIVE_ASSEMBLY}} {params.out_prefix}  1>{log.std} 2>&1;"


