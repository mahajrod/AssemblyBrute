localrules: gather_pretext_track_info

def select_tracks(wildcards):

    stage_parameter_dict = parse_stage_parameters_from_path(wildcards.fasta_dir)
    stage = stage_parameter_dict["stage"]
    parameter_label = stage_parameter_dict["parameters"]
    haplotype_list = stage_dict[stage].parameters[parameter_label]["haplotype_list"]

    track_prefix = "{0}/assembly_qc/tracks/{1}/{1}".format(wildcards.fasta_dir,
                                                           wildcards.fasta_prefix)
    track_dict = {"gap": f"{track_prefix}.gap.track.bedgraph", }
    if not config["skip_telomere"]:
        track_ext = parameters["tool_options"]["pretextview"]["track_ext"]["telomere"]
        if not config["skip_telomere_container"]:
            track_dict[f"cannon_telomere.{track_ext}"] = f"{track_prefix}.canonical.telomere.pretext.bedgraph"
            track_dict[f"non_cannon_telomere.{track_ext}"] = f"{track_prefix}.non_canonical.telomere.pretext.bedgraph"
        if not config["skip_telomere_tidk"]:
            track_dict[f"canon_telomere_tidk.{track_ext}"] = f"{track_prefix}.canonical_tidk_telomere_filtered.win%s.step%s.track.bedgraph" % (parameters["tool_options"]["assembly_qc"]["telomere_tidk_search"]["window_size"],
                                                                                                                                                        parameters["tool_options"]["assembly_qc"]["telomere_tidk_search"]["window_size"])
            track_dict[f"non_canon_telomere_tidk.{track_ext}"] = f"{track_prefix}.non_canonical_tidk_telomere_filtered.win%s.step%s.track.bedgraph" % (parameters["tool_options"]["assembly_qc"]["telomere_tidk_search"]["window_size"],
                                                                                                                                                                parameters["tool_options"]["assembly_qc"]["telomere_tidk_search"]["window_size"])

    window_track_suffix_dict = {}
    if not config["skip_pretext_10k_1k_tracks"]:
        window_track_suffix_dict["10k_1k"] = "win10000.step1000"
    if not config["skip_pretext_100k_10k_tracks"]:
        window_track_suffix_dict["100k_10k"] = "win100000.step10000"
    if not config["skip_pretext_1000k_100k_tracks"]:
        window_track_suffix_dict["1000k_100k"] = "win1000000.step100000"

    track_type_dict = {"gc":           "gc",
                       "windowmasker": "windowmasker"}
    track_ext_dict = {"gc":              parameters["tool_options"]["pretextview"]["track_ext"]["gc"],
                      "windowmasker":    parameters["tool_options"]["pretextview"]["track_ext"]["windowmasker"]}

    if not config["skip_trf"]:
        track_type_dict["trf"] = "trf"
        track_ext_dict["trf"] = parameters["tool_options"]["pretextview"]["track_ext"]["trf"]
    if not config["skip_pretext_coverage_tracks"]:
        for datatype in config["data_feature_dict"]["pretext_coverage_track"]:
            if datatype in config["data_feature_dict"]["pretext_per_hap_track"]:
                #Inject per haplotype tracks, i.e. tracks calculated for each haplotype individually and later labeled for merged haplotype
                for haplotype in haplotype_list:
                    track_type_dict[f"{haplotype}@{datatype}_all_nodup_mean"] = f"{haplotype}@{datatype}_all_nodup_reads_mean_coverage"
                    track_ext_dict[f"{haplotype}@{datatype}_all_nodup_mean"] = parameters["tool_options"]["pretextview"]["track_ext"]["coverage"]
            else:
                track_type_dict[f"{datatype}_all_nodup_mean"] = f"{datatype}_all_nodup_reads_mean_coverage"
                track_ext_dict[f"{datatype}_all_nodup_mean"] = parameters["tool_options"]["pretextview"]["track_ext"]["coverage"]
        for datatype in config["ext_data_feature_dict"]:
            for track_name in config["ext_data_feature_dict"][datatype]["pretext_coverage_track"]:
                if track_name in config["ext_data_feature_dict"][datatype]["pretext_per_hap_track"]:
                    #Inject per haplotype tracks, i.e. tracks calculated for each haplotype individually and later labeled for merged haplotype
                    for haplotype in haplotype_list:
                        track_type_dict[f"{haplotype}@ext@{datatype}@{track_name}_all_nodup_mean"] = f"{haplotype}@ext@{datatype}@{track_name}_all_nodup_reads_mean_coverage"
                        track_ext_dict[f"{haplotype}@ext@{datatype}@{track_name}_all_nodup_mean"] = parameters["tool_options"]["pretextview"]["track_ext"]["coverage"]

                        track_type_dict[f"{haplotype}@ext@{datatype}@{track_name}_hq_mapping_mean"] = f"{haplotype}@ext@{datatype}@{track_name}_hq_mapping_mean_coverage"
                        track_ext_dict[f"{haplotype}@ext@{datatype}@{track_name}_hq_mapping_mean"] = parameters["tool_options"]["pretextview"]["track_ext"]["coverage"]
                else:
                    track_type_dict[f"ext@{datatype}@{track_name}_all_nodup_mean"] = f"ext@{datatype}@{track_name}_all_nodup_reads_mean_coverage"
                    track_ext_dict[f"ext@{datatype}@{track_name}_all_nodup_mean"] = parameters["tool_options"]["pretextview"]["track_ext"]["coverage"]

    for suffix in window_track_suffix_dict:
        for track_type in track_type_dict:
            track_dict[f"{track_type}_{suffix}.{track_ext_dict[track_type]}"] = f"{track_prefix}.{track_type_dict[track_type]}.{window_track_suffix_dict[suffix]}.track.bedgraph"

    if not config["skip_busco"]:
        track_ext = parameters["tool_options"]["pretextview"]["track_ext"]["busco"]
        for busco_lineage in config["busco_lineage_list"]:
            for busco_type in ["single_copy", "duplicated", "fragmented"]:
                track_dict[f"BUSCO.{busco_lineage}.{busco_type}.{track_ext}"] = f"{track_prefix}.{busco_lineage}.busco5.{busco_type}.track.bedgraph"
    if not config["skip_purge_dups_qc"]:
        track_ext = parameters["tool_options"]["pretextview"]["track_ext"]["purge_dups"]
        for datatype in parameters["tool_options"]["assembly_qc"]["purge_dups"]["datatype_list"]:
            for artefact_type in ["junk", "ovlp", "haplotig", "repeat", "highcov"]:
                track_dict[f"purge_dups.{datatype}.{artefact_type}.{track_ext}"] = f"{track_prefix}.purge_dups.{datatype}.{artefact_type}.track.bedgraph"
    #print(track_dict)
    return track_dict

rule gather_pretext_track_info:
    input:
        unpack(select_tracks),
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        track_info="{fasta_dir}/{fasta_prefix}.pretext.track.info",
    log:
        log="{fasta_dir}/log/gather_pretext_track_info.{fasta_prefix}.log",
    run:
        track_dict = dict(input.items())

        track_dict.pop("log_dir")
        with open(log.log, "w") as log_fd:
            with redirect_stdout(log_fd), redirect_stderr(log_fd):
                with open(output.track_info, "w") as track_info_fd:
                    for track in track_dict:
                        track_info_fd.write(f"{track}\t{track_dict[track]}\n")

rule pretext_inject_tracks:
    input:
        map="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.{subset}.rmdup.mapq{mapq}.{pretext_res}.pretext",
        filtered_out="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.{subset}.rmdup.mapq{mapq}.{pretext_res}.filtered_out.ids",
        track_info="{fasta_dir}/{fasta_prefix}.pretext.track.info",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        updated_map="{fasta_dir}/{fasta_prefix}/alignment/{phasing_kmer_length}/{fasta_prefix}.{phasing_kmer_length}.{subset, all}.rmdup.mapq{mapq}.{pretext_res}.tracks.pretext",
    log:
        injection="{fasta_dir}/log/pretext_inject_tracks.{fasta_prefix}.{phasing_kmer_length}.{subset}.{mapq}.{pretext_res}.injection.log",
        cluster_log="{fasta_dir}/log/pretext_inject_tracks.{fasta_prefix}.{phasing_kmer_length}.{subset}.{mapq}.{pretext_res}.cluster.log",
        cluster_err="{fasta_dir}/log/pretext_inject_tracks.{fasta_prefix}.{phasing_kmer_length}.{subset}.{mapq}.{pretext_res}.cluster.err"
    benchmark:
        "{fasta_dir}/log/pretext_inject_tracks.{fasta_prefix}.{phasing_kmer_length}.{subset}.{mapq}.{pretext_res}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("pretext_inject_tracks"),
        cpus=parameters["threads"]["pretext_inject_tracks"] ,
        time=parameters["time"]["pretext_inject_tracks"],
        mem=parameters["memory_mb"]["pretext_inject_tracks"]
    threads: parameters["threads"]["pretext_inject_tracks"]

    shell:
        " cp -f {input.map} {output.updated_map}; "
        " > {log.injection}; "
        " TRACK_NAME_ARRAY=(`cut -f 1 {input.track_info}`); "
        " TRACK_PATH_ARRAY=(`cut -f 2 {input.track_info}`); "
        " TRACK_NUMBER=${{#TRACK_NAME_ARRAY[@]}}; "
        " for ((INDEX=0; INDEX<${{#TRACK_NAME_ARRAY[@]}}; INDEX++)); "
        " do "
        "     FILTERED_TRACK=${{TRACK_PATH_ARRAY[INDEX]}}.{wildcards.phasing_kmer_length}.{wildcards.subset}.{wildcards.mapq}.{wildcards.pretext_res}.tmp; "
        "     workflow/scripts/curation/filter_bed_by_scaffolds.py -d {input.filtered_out} "
        "         -i ${{TRACK_PATH_ARRAY[INDEX]}} > ${{FILTERED_TRACK}} 2>>{log.injection}; "
        "     if [ -s ${{FILTERED_TRACK}} ]; "
        "     then "
        "         echo \"Injecting track ${{TRACK_NAME_ARRAY[INDEX]}}: ${{FILTERED_TRACK}}\" >> {log.injection}; "
        "         cat ${{FILTERED_TRACK}} | awk '{{printf \"%s\\t%i\\t%i\\t%i\\n\",$1,$2,$3,$4}}' | "
        "             PretextGraph -i {output.updated_map} -n ${{TRACK_NAME_ARRAY[INDEX]}} >> {log.injection} 2>&1;"
        "     fi; "
        "     rm ${{FILTERED_TRACK}}; " 
        " done "
