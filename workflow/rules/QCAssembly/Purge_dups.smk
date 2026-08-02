
ruleorder: create_purge_dups_track_for_combined_haplotype > copy_purge_dups_track

"""
rule get_purged_seqs: #
    input:
        raw_dups_bed=config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}/purge_dups/{datatype}/{genome_prefix}.{assembly_stage}.{haplotype}.dups.raw.bed",
        assembly=config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta",
        log_dir=ancient(config["out_dir"] / "{assembly_stage}/{parameters}/log/")
    output:
        filtered_bed=config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}/purge_dups/{datatype}/{genome_prefix}.{assembly_stage}.{haplotype}.dups.filtered.bed",
        purged=config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}/purge_dups/{datatype}/{genome_prefix}.{assembly_stage}.{haplotype}.purged.fasta",
        hapdups=config["out_dir"] / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}/purge_dups/{datatype}/{genome_prefix}.{assembly_stage}.{assembly_stage}.{haplotype}.hap.fasta",
    params:
        blacklist_option=lambda wildcards: parse_option("purging_blacklist",
                                                        stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["option_set"],
                                                        option_prefix="-b", expression=lambda l: ",".join(l)),
        whitelist_option=lambda wildcards: parse_option("purging_whitelist",
                                                        stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["option_set"],
                                                        option_prefix="-w", expression=lambda l: ",".join(l)),
    log:
        get_seqs=config["out_dir"] / "{assembly_stage}/{parameters}/log/qc_get_purged_seqs.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.get_seqs.log",
        filter=config["out_dir"] / "{assembly_stage}/{parameters}/log/qc_get_purged_seqs.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.filter.log",
        ln=config["out_dir"] / "{assembly_stage}/{parameters}/log/qc_get_purged_seqs.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.ln.log",
        cluster_log=config["out_dir"] / "{assembly_stage}/{parameters}/log/qc_get_purged_seqs.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.cluster.log",
        cluster_err=config["out_dir"] / "{assembly_stage}/{parameters}/log/qc_get_purged_seqs.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.cluster.err"
    benchmark:
        config["out_dir"] / "{assembly_stage}/{parameters}/log/qc_get_purged_seqs.{parameters}.{genome_prefix}.{assembly_stage}.{haplotype}.{datatype}.benchmark.txt"
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
        " workflow/scripts/purge_dups/filter_dups_bed.py -i {input.raw_dups_bed} "
        " {params.blacklist_option} {params.whitelist_option} -o {output.filtered_bed} > {log.filter} 2>&1; "
        " OUT_DIR=`dirname {output.filtered_bed}`; "
        " PURGE_DUPS_BED=`realpath -s {output.filtered_bed}`; "
        " REFERENCE=`realpath -s {input.assembly}`; "
        " GET_SEQ_LOG=`realpath -s {log.get_seqs}`; "
        " LN_LOG=`realpath -s {log.ln}`; "
        " cd ${{OUT_DIR}}; "
        " get_seqs -p {wildcards.genome_prefix}.{wildcards.assembly_stage}.{wildcards.haplotype} ${{PURGE_DUPS_BED}} ${{REFERENCE}} > ${{GET_SEQ_LOG}} 2>&1; "
        " for FILE in *.fa; do mv ${{FILE}} ${{FILE%fa}}fasta; done; "
"""

rule copy_purge_dups_track:
    priority: 500
    input:
        bedgraph="{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/{fasta_prefix}.purge_dups.{datatype}.{artefact_type}.track.bedgraph",
        bed="{fasta_dir}/{fasta_prefix}/purge_dups/{datatype}/{fasta_prefix}.purge_dups.{datatype}.{artefact_type}.track.bed",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.purge_dups.{datatype}.{artefact_type}.track.bedgraph",
        bed="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.purge_dups.{datatype}.{artefact_type}.track.bed",
    log:
        log="{fasta_dir}/log/create_purge_dups_track.{fasta_prefix}.{datatype}.{artefact_type}.single_copy.log",
        cluster_log="{fasta_dir}/log/create_purge_dups_track.{fasta_prefix}.{datatype}.{artefact_type}.cluster.log",
        cluster_err="{fasta_dir}/log/create_purge_dups_track.{fasta_prefix}.{datatype}.{artefact_type}.cluster.err"
    benchmark:
        "{fasta_dir}/log/create_purge_dups_track.{fasta_prefix}.{datatype}.{artefact_type}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("busco5_intersect_all"),
        cpus=parameters["threads"]["busco5_intersect_all"],
        time=parameters["time"]["busco5_intersect_all"],
        mem=parameters["memory_mb"]["busco5_intersect_all"],
    threads:
        parameters["threads"]["busco5_intersect_all"]
    shell:
        " cp -f  {input.bed} {input.bedgraph} `dirname {output.bedgraph}` 2>{log.log}; "


rule create_purge_dups_track_for_combined_haplotype:
    priority: 500
    input:
        single_haplotype_tracks=lambda wildcards: expand(config["out_dir"] / ("%s/%s/assembly_qc/tracks/%s.%s.{haplotype}/%s.%s.{haplotype}.purge_dups.%s.%s.track.bed" % (wildcards.assembly_stage,
                                                                                                                                                                           wildcards.parameters,
                                                                                                                                                                           wildcards.genome_prefix,
                                                                                                                                                                           wildcards.assembly_stage,
                                                                                                                                                                           wildcards.genome_prefix,
                                                                                                                                                                           wildcards.assembly_stage,
                                                                                                                                                                           wildcards.datatype,
                                                                                                                                                                           wildcards.artefact_type)),
                                                           haplotype=stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"],
                                                           allow_missing=True),
        log_dir=ancient(config["out_dir"] / "{assembly_stage}/{parameters}/log/"),
    output:
        combined_track=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{merged_haplotype}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.purge_dups.{datatype}.{artefact_type}.track.bedgraph",
    params:
        haplotype_list=lambda wildcards: stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"],
        dir_prefix=lambda wildcards: config["out_dir"] / ("%s/%s/assembly_qc/tracks/%s.%s" % (wildcards.assembly_stage,
                                                                                         wildcards.parameters,
                                                                                         wildcards.genome_prefix,
                                                                                         wildcards.assembly_stage)),
        track_prefix=lambda wildcards: "%s.%s" % (wildcards.genome_prefix, wildcards.assembly_stage),
        suffix=lambda wildcards: "purge_dups.%s.%s" % (wildcards.datatype, wildcards.artefact_type),
    log:
        log=config["out_dir"] / "{assembly_stage}/{parameters}/log/create_purge_dups_track_for_combined_haplotype.{assembly_stage}.{parameters}.{genome_prefix}.{merged_haplotype}.purge_dups.{datatype}.{artefact_type}.log",
        cluster_log=config["out_dir"] / "{assembly_stage}/{parameters}/log/create_purge_dups_track_for_combined_haplotype.{assembly_stage}.{parameters}.{genome_prefix}.{merged_haplotype}.purge_dups.{datatype}.{artefact_type}.cluster.log",
        cluster_err=config["out_dir"] / "{assembly_stage}/{parameters}/log/create_purge_dups_for_combined_haplotype.{assembly_stage}.{parameters}.{genome_prefix}.{merged_haplotype}.purge_dups.{datatype}.{artefact_type}.cluster.err"
    benchmark:
        config["out_dir"] / "{assembly_stage}/{parameters}/log/create_purge_dups_for_combined_haplotype.{assembly_stage}.{parameters}.{genome_prefix}.{merged_haplotype}.purge_dups.{datatype}.{artefact_type}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("busco5_intersect_all"),
        cpus=parameters["threads"]["busco5_intersect_all"],
        time=parameters["time"]["busco5_intersect_all"],
        mem=parameters["memory_mb"]["busco5_intersect_all"],
    threads:
        parameters["threads"]["busco5_intersect_all"]
    shell:
        " > {log.log}; "
        " > {output.combined_track}; "
        " for HAP in {params.haplotype_list}; "
        " do "
        "     echo \"Processing ${{HAP}}...\" >> {log.log} 2>&1; "
        "     HAP_TRACK={params.dir_prefix}.${{HAP}}/{params.track_prefix}.${{HAP}}.{params.suffix}.track.bedgraph; "
        "     sed 's/^/'${{HAP}}'\./' ${{HAP_TRACK}} >> {output.combined_track} 2>>{log.log}; "
        " done;"
