localrules: create_final_combined_dedup_links

rule combined_dedup: #rule for combined deduplication using both purge_dups and hapsolo. Only contigs marked by both tools will be removed
    input:
        raw_purge_dups_bed=config["out_dir"] / "dedup/{parameters_prefix}@{busco_lineage}/{input_fasta_prefix}/{input_fasta_prefix}/purge_dups/{datatype}/{input_fasta_prefix}.dups.raw.bed",
        hapsolo_duplication_ids=config["out_dir"] / "dedup/{parameters_prefix}@{busco_lineage}/{input_fasta_prefix}/{input_fasta_prefix}/hapsolo_qc/{busco_lineage}/{input_fasta_prefix}.dups.ids",
        reference=config["out_dir"] / "dedup/{parameters_prefix}@{busco_lineage}/{input_fasta_prefix}/{input_fasta_prefix}.fasta",
        log_dir=ancient(config["out_dir"] / "dedup/{parameters_prefix}@{busco_lineage}/log/"),
    output:
        final_dups_filtered_bed=config["out_dir"] / "dedup/{parameters_prefix, [^/]*combo_purge[^/]*}@{busco_lineage}/{input_fasta_prefix, [^/]*input[^/]*}/{input_fasta_prefix}/combo_purge/{datatype}.{busco_lineage}/{input_fasta_prefix}.dups.final.bed",
        purged=config["out_dir"] / "dedup/{parameters_prefix, [^/]*combo_purge[^/]*}@{busco_lineage}/{input_fasta_prefix, [^/]*input[^/]*}/{input_fasta_prefix}/combo_purge/{datatype}.{busco_lineage}/{input_fasta_prefix}.purged.fasta",
        hapdups=config["out_dir"] / "dedup/{parameters_prefix, [^/]*combo_purge[^/]*}@{busco_lineage}/{input_fasta_prefix, [^/]*input[^/]*}/{input_fasta_prefix}/combo_purge/{datatype}.{busco_lineage}/{input_fasta_prefix}.hap.fasta",
    params:
        blacklist_option=lambda wildcards: parse_option("purging_blacklist",
                                                        stage_dict["dedup"].parameters[wildcards.parameters_prefix + "@" + wildcards.busco_lineage]["option_set"],
                                                        option_prefix="-b", expression=lambda l: ",".join(l)),
        whitelist_option=lambda wildcards: parse_option("purging_whitelist",
                                                        stage_dict["dedup"].parameters[wildcards.parameters_prefix + "@" + wildcards.busco_lineage]["option_set"],
                                                        option_prefix="-w", expression=lambda l: ",".join(l)),
        rel_ref_path=lambda wildcards: get_relative_path(config["out_dir"] / "dedup/{0}@{2}/{1}/{1}.fasta".format(wildcards.parameters_prefix,
                                                                                                                  wildcards.input_fasta_prefix,
                                                                                                                  wildcards.busco_lineage),
                                                         config["out_dir"] / "dedup/{0}@{3}/{1}/{1}/combo_purge/{2}.{3}/{1}.purged.fasta".format(wildcards.parameters_prefix,
                                                                                                                                                 wildcards.input_fasta_prefix,
                                                                                                                                                 wildcards.datatype,
                                                                                                                                                 wildcards.busco_lineage))
    log:
        get_seqs=config["out_dir"] / "dedup/{parameters_prefix}@{busco_lineage}/log/combined_dedup.{parameters_prefix}@{busco_lineage}.{input_fasta_prefix}.{datatype}.get_seqs.log",
        filter=config["out_dir"] / "dedup/{parameters_prefix}@{busco_lineage}/log/combined_dedup.{parameters_prefix}@{busco_lineage}.{input_fasta_prefix}.{datatype}.filter.log",
        ln=config["out_dir"] / "dedup/{parameters_prefix}@{busco_lineage}/log/combined_dedup.{parameters_prefix}@{busco_lineage}.{input_fasta_prefix}.{datatype}.ln.log",
        cluster_log=config["out_dir"] / "dedup/{parameters_prefix}@{busco_lineage}/log/combined_dedup.{parameters_prefix}@{busco_lineage}.{input_fasta_prefix}.{datatype}.cluster.log",
        cluster_err=config["out_dir"] / "dedup/{parameters_prefix}@{busco_lineage}/log/combined_dedup.{parameters_prefix}@{busco_lineage}.{input_fasta_prefix}.{datatype}.cluster.err"
    benchmark:
        config["out_dir"] / "dedup/{parameters_prefix}@{busco_lineage}/log/combined_dedup.{parameters_prefix}@{busco_lineage}.{input_fasta_prefix}.{datatype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("purge_dups"),
        cpus=parameters["threads"]["purge_dups"] ,
        time=parameters["time"]["purge_dups"],
        mem=parameters["memory_mb"]["purge_dups"]
    threads: parameters["threads"]["purge_dups"]

    shell:
        " workflow/scripts/purge_dups/filter_dups_bed.py -i {input.raw_purge_dups_bed} {params.blacklist_option} "
        "         {params.whitelist_option} --scaffold_whitelist {input.hapsolo_duplication_ids} "
        "         -o {output.final_dups_filtered_bed} > {log.filter} 2>&1; "
        " OUT_DIR=`dirname {output.final_dups_filtered_bed}`; "
        " PURGE_DUPS_BED=`realpath -s {output.final_dups_filtered_bed}`; "
        " REFERENCE=`realpath -s {input.reference}`; "
        " GET_SEQ_LOG=`realpath -s {log.get_seqs}`; "
        " LN_LOG=`realpath -s {log.ln}`; "
        " cd ${{OUT_DIR}}; "
        " if [ -s ${{PURGE_DUPS_BED}}  ];" # check if ${{PURGE_DUPS_BED}} exists and is not empty
        " then "
        "     get_seqs -p {wildcards.input_fasta_prefix} ${{PURGE_DUPS_BED}} ${{REFERENCE}} > ${{GET_SEQ_LOG}} 2>&1; "
        "     for FILE in *.fa; do mv ${{FILE}} ${{FILE%fa}}fasta; done; "
        " else "
        "     cp -f {params.rel_ref_path} `basename {output.purged}` >> ${{GET_SEQ_LOG}} 2>&1; "
        "     touch {wildcards.input_fasta_prefix}.hap.fasta; "
        " fi; "


use rule create_local_links as create_final_combined_dedup_links with:
    input:
        fasta=lambda wildcards: config["out_dir"] / ("dedup/{parameters}/{genome_prefix}.input.{haplotype}/{genome_prefix}.input.{haplotype}/combo_purge/%s.%s/{genome_prefix}.input.{haplotype}.purged.fasta" % ("_".join(stage_dict["dedup"].parameters[wildcards.parameters]["option_set"]["main_datatypes"]),
                                                                                                                                                                                                                     wildcards.parameters.split("@")[-1])),

        log_dir=ancient(config["out_dir"] / "dedup/{parameters}/log/")
    output:
        fasta=config["out_dir"] / "dedup/{parameters, [^/]*combo_purge[^/]*}/{genome_prefix}.dedup.{haplotype, hap[^_./@]+}.fasta",
    log:
        ln=config["out_dir"] / "dedup/{parameters}/log/create_final_combined_dedup_links.{parameters}.{genome_prefix}.dedup.{haplotype}.ln.log",
