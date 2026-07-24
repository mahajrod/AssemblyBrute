
def get_meryl_db_for_merqury(wildcards):
    datatype_list = wildcards.datatype.split("_")
    kmer_datatype_list = []
    for datatype in datatype_list:
        if datatype in config["data"]:
            kmer_datatype_list.append(datatype)
        else:
            print(f"WARNING. QC datatype {datatype} is missing in the input. Excluding it from the QC!!!")


    return config["out_dir"] / "kmer/{0}/final/{0}.final.{1}.meryl".format("_".join(kmer_datatype_list),
                                                                               config["final_kmer_length"],)

rule merqury:
    input:
        meryl_db_dir = get_meryl_db_for_merqury,
        primary_assembly=lambda wildcards: config["out_dir"] / "{0}/{1}/{2}.{0}.{3}.fasta".format(wildcards.assembly_stage,
                                                                                             wildcards.parameters,
                                                                                             wildcards.genome_prefix,
                                                                                             stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"][0]),
        alternative_assembly=lambda wildcards: config["out_dir"] / "{0}/{1}/{2}.{0}.{3}.fasta".format(wildcards.assembly_stage,
                                                                                             wildcards.parameters,
                                                                                             wildcards.genome_prefix,
                                                                                             "hap2") if stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"][0] != "hap0" else [],
    output:
        qv_file=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/merqury/{datatype}/{genome_prefix}.{assembly_stage}.{datatype}.qv",
        completeness_stats_file=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/merqury/{datatype}/{genome_prefix}.{assembly_stage}.{datatype}.completeness.stats",
    params:
        out_prefix=lambda wildcards: "{0}.{1}.{2}".format(wildcards.genome_prefix,
                                                          wildcards.assembly_stage,
                                                          wildcards.datatype)
    log:
        std=config["out_dir"].resolve() / "log/merqury.{assembly_stage}.{parameters}.{genome_prefix}.{datatype}.log",
        mkdir_log=(config["out_dir"]).resolve() / "log/merqury.{assembly_stage}.{parameters}.{genome_prefix}.{datatype}.mkdir.log",
        cd_log=(config["out_dir"]).resolve() / "log/merqury.{assembly_stage}.{parameters}.{genome_prefix}.{datatype}.cd.log",
        cluster_log=(config["out_dir"]).resolve() / "log/merqury.{assembly_stage}.{parameters}.{genome_prefix}.{datatype}.cluster.log",
        cluster_err=(config["out_dir"]).resolve() / "log/merqury.{assembly_stage}.{parameters}.{genome_prefix}.{datatype}.cluster.err"
    benchmark:
        config["out_dir"] / "log/merqury.{assembly_stage}.{parameters}.{genome_prefix}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["kmer"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["kmer"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
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
         " OMP_NUM_THREADS={threads} merqury.sh ${{MERYL_DB}} ${{PRIMARY_ASSEMBLY}} ${{ALTERNATIVE_ASSEMBLY}} {params.out_prefix} >{log.std} 2>&1;"
