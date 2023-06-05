rule merqury: # TODO: add handling for cases of haploid and polyploid genomes
    input:
        meryl_db_dir=output_dict["kmer"] / "{0}/filtered/{0}.filtered.{1}.meryl".format(config["final_kmer_datatype"],
                                                                                        config["final_kmer_length"],),
        primary_assembly=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.hap1.fasta",
        alternative_assembly=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.hap2.fasta",
    output:
        #dir=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/merqury/",
        qv_file=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/merqury/{genome_prefix}.{assembly_stage}.qv",
        completeness_stats_file=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/merqury/{genome_prefix}.{assembly_stage}.completeness.stats",
    params:
        dir=lambda wildcards: out_dir_path / "{0}/{1}/assembly_qc/merqury/".format(wildcards.assembly_stage,
                                                                                   wildcards.parameters),
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
         " MERYL_DB=`realpath {input.meryl_db_dir}`;"
         " PRIMARY_ASSEMBLY=`realpath -s {input.primary_assembly}`;"
         " ALTERNATIVE_ASSEMBLY=`realpath -s {input.alternative_assembly}`;"
         " cd {params.dir}; "
         " OMP_NUM_THREADS={threads} merqury.sh ${{MERYL_DB}} "
         " ${{PRIMARY_ASSEMBLY}} ${{ALTERNATIVE_ASSEMBLY}} {params.out_prefix}  1>{log.std} 2>&1;"
         #" OMP_NUM_THREADS={threads} merqury.sh {input.meryl_db_dir} "
         #" {input.primary_assembly} {input.alternative_assembly} {params.out_prefix}  1>{log.std} 2>&1 || true;"


