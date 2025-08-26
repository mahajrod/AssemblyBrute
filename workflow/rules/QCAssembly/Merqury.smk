
#stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"]

def get_meryl_db_for_merqury(wildcards):
    if ("qc_datatypes" in config) and config["qc_datatypes"]:
        qc_datatypes = config["qc_datatypes"]
    else:
        qc_datatypes = stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["option_set"]["qc_datatypes"]
    print("QC datatypes: {0}\n".format('.'.join(qc_datatypes)))
    kmer_datatype_list = []
    filtered_flag = False
    for datatype in qc_datatypes:
        if datatype in input_filedict:
            kmer_datatype_list.append(datatype)
            if datatype in config["filtered_data"]:
                filtered_flag = True

    #print(output_dict["kmer"] / "{0}/{1}/{0}.{1}.{2}.meryl".format("_".join(kmer_datatype_list),
    #                                                                "filtered" if filtered_flag else "raw",
    #                                                                 config["final_kmer_length"],))
    return output_dict["kmer"] / "{0}/{1}/{0}.{1}.{2}.meryl".format("_".join(kmer_datatype_list),
                                                                    "filtered" if filtered_flag else "raw",
                                                                     config["final_kmer_length"],)


rule merqury: # TODO: add handling for cases of haploid and polyploid genomes
    input:
        #meryl_db_dir=output_dict["kmer"] / "{0}/{1}/{0}.{1}.{2}.meryl".format(config["final_kmer_datatype"],
        #                                                                      "filtered" if config["final_kmer_datatype"] in config["filtered_data"] else "raw",
        #                                                                      config["final_kmer_length"],) ,
        meryl_db_dir = get_meryl_db_for_merqury,
        primary_assembly=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(wildcards.assembly_stage,
                                                                                             wildcards.parameters,
                                                                                             wildcards.genome_prefix,
                                                                                             stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["haplotype_list"][0]),
        alternative_assembly=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(wildcards.assembly_stage,
                                                                                             wildcards.parameters,
                                                                                             wildcards.genome_prefix,
                                                                                             "hap2") if stage_dict[wildcards.assembly_stage]["parameters"][wildcards.parameters]["haplotype_list"][0] != "hap0" else [],
    output:
        qv_file=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/merqury/{genome_prefix, [^/]+}.{assembly_stage}.qv",
        completeness_stats_file=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/merqury/{genome_prefix, [^/]+}.{assembly_stage}.completeness.stats",
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
         " MERQURY_SCRIPT=`realpath workflow/external_tools/merqury/merqury.sh`; "
         " MERQURY=`dirname ${{MERQURY_SCRIPT}}`; "
         " MERQURY_DIR=`dirname ${{MERQURY_SCRIPT}}`; "
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
         " PATH=${{MERQURY_DIR}}:${{PATH}} OMP_NUM_THREADS={threads} ${{MERQURY_SCRIPT}} ${{MERYL_DB}} "
         " ${{PRIMARY_ASSEMBLY}} ${{ALTERNATIVE_ASSEMBLY}} {params.out_prefix}  1>{log.std} 2>&1;"


