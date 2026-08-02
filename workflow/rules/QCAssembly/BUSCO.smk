localrules: busco_intersect_haplotypes, busco_intersect_stages
localrules: create_busco_tracks, create_busco_tracks_for_combined_haplotype
ruleorder: create_busco_tracks_for_combined_haplotype > create_busco_tracks

rule busco_intersect_haplotypes: # Downloading of busco datasets is performed by a different rule to avoid conflict between different instances of busco5
    priority: 500
    input:
        busco_tables=lambda wildcards: expand(config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}.{busco_version}.full_table.tsv",
                                              haplotype=stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"], # if wildcards.assembly_stage != "draft_qc" else haplotype_list,
                                              allow_missing=True,)
    params:
        haplotypes=lambda wildcards: ",".join(stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"])
    output:
        busco_legend=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.legend",
        busco_orderlist=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.orderlist",
        busco_merged_tsv=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.merged.tsv",
        busco_len=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.len",
        busco_counts_bedgraph=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.counts.bedgraph",
        busco_counts_png=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.png",
        no_complete_busco_len=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.no_complete.len",
        no_complete_busco_counts_bedgraph=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.no_complete.counts.bedgraph",
        no_complete_busco_counts_png=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.no_complete.png",
        informative_busco_len=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.informative.len",
        informative_busco_counts_bedgraph=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.informative.counts.bedgraph",
        informative_busco_counts_png=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.informative.png",
    log:
        std=config["out_dir"] / "log/busco_intersect_haplotypes.{assembly_stage}.{parameters}.{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.log",
        draw=config["out_dir"] / "log/busco_intersect_haplotypes.{assembly_stage}.{parameters}.{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.draw.log",
        draw_no_complete=config["out_dir"] / "log/busco_intersect_haplotypes.{assembly_stage}.{parameters}.{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.draw_no_complete.log",
        draw_informative=config["out_dir"] / "log/busco_intersect_haplotypes.{assembly_stage}.{parameters}.{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.draw_informative.log",
        cluster_log=config["out_dir"] / "log/busco_intersect_haplotypes.{assembly_stage}.{parameters}.{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.cluster.log",
        cluster_err=config["out_dir"] / "log/busco_intersect_haplotypes.{assembly_stage}.{parameters}.{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.cluster.err"
    benchmark:
        config["out_dir"] / "log/busco5_intersect_haplotypes.{assembly_stage}.{parameters}.{genome_prefix}.{assembly_stage}.{busco_lineage}.{busco_version}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("busco5_intersect_haplotypes"),
        cpus=parameters["threads"]["busco5_intersect_haplotypes"],
        time=parameters["time"]["busco5_intersect_haplotypes"],
        mem=parameters["memory_mb"]["busco5_intersect_haplotypes"],
    threads:
        parameters["threads"]["busco5_intersect_haplotypes"]
    shell:
         " OUTPUT_PREFIX={output.busco_legend};"
         " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.legend}};  "
         " workflow/scripts/busco/intersect_busco_results.py -b `echo '{input.busco_tables}' | tr ' ' ',' ` "
         "     -l {params.haplotypes} -o ${{OUTPUT_PREFIX}} > {log.std} 2>&1; "
         " draw_features.py  -i {output.busco_counts_bedgraph} -o ${{OUTPUT_PREFIX}} -t bedgraph  "
         "     -n {output.busco_len} -z {output.busco_orderlist} -g {output.busco_legend} --hide_track_label "
         "     --color_column_name value -l {wildcards.busco_lineage} --figure_header_height 2 --subplots_adjust_top 0.7 "
         "     --x_tick_type int_number > {log.draw} 2>&1; "
         " draw_features.py  -i {output.no_complete_busco_counts_bedgraph} -o ${{OUTPUT_PREFIX}}.no_complete "
         "     -t bedgraph  -n {output.no_complete_busco_len} -z {output.busco_orderlist} -g {output.busco_legend} "
         "     --hide_track_label  --color_column_name value -l {wildcards.busco_lineage} "
         "     --figure_header_height 2 --subplots_adjust_top 0.7 "
         "     --x_tick_type int_number > {log.draw_no_complete} 2>&1; "
         " draw_features.py  -i {output.informative_busco_counts_bedgraph} -o ${{OUTPUT_PREFIX}}.informative "
         "     -t bedgraph  -n {output.informative_busco_len} -z {output.busco_orderlist} -g {output.busco_legend} "
         "     --hide_track_label  --color_column_name value -l {wildcards.busco_lineage} "
         "     --figure_header_height 2 --subplots_adjust_top 0.7 "
         "     --x_tick_type int_number > {log.draw_informative} 2>&1; "
         " touch {output.busco_counts_png}; "
         " touch {output.no_complete_busco_counts_png}; "
         " touch {output.informative_busco_counts_png}"


def get_busco_table_for_all_assemblies_in_chain_per_haplotype(wildcards):
    busco_table_list = []
    parameters_dict = get_parameters_for_all_stages_in_chain(wildcards.parameters)
    for skip_stage in config["skip_qc_fasta_stage_list"]:
        if skip_stage in parameters_dict:
            parameters_dict.pop(skip_stage)
    for stage in parameters_dict:
        busco_table_list += expand(config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}.{busco_version}.full_table.tsv",
                                               assembly_stage=[stage],
                                               parameters=[parameters_dict[stage]],
                                               allow_missing=True,)
    return  busco_table_list

rule busco_intersect_stages:
    priority: 500
    input:
        busco_tables=get_busco_table_for_all_assemblies_in_chain_per_haplotype,
    params:
        stages=lambda wildcards: ",".join(set(get_parameters_for_all_stages_in_chain(wildcards.parameters).keys()) - set(config["skip_qc_fasta_stage_list"])),
    output:
        busco_legend=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.legend",
        busco_orderlist=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.orderlist",
        busco_merged_tsv=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.merged.tsv",
        busco_informative_tsv=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.informative.tsv",
        busco_len=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.len",
        busco_counts_bedgraph=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.counts.bedgraph",
        busco_counts_png=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.png",
        no_complete_busco_len=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.no_complete.len",
        no_complete_busco_counts_bedgraph=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.no_complete.counts.bedgraph",
        no_complete_busco_counts_png=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.no_complete.png",
        informative_busco_len=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.informative.len",
        informative_busco_counts_bedgraph=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.informative.counts.bedgraph",
        informative_busco_counts_png=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.informative.png",
    log:
        std=config["out_dir"] / "log/busco_intersect_stages.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.log",
        draw=config["out_dir"] / "log/busco_intersect_stages.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.draw.log",
        draw_no_complete=config["out_dir"] / "log/busco_intersect_stages.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.draw_no_complete.log",
        draw_informative=config["out_dir"] / "log/busco_intersect_stages.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.draw_informative.log",
        cluster_log=config["out_dir"] / "log/busco_intersect_stages.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.cluster.log",
        cluster_err=config["out_dir"] / "log/busco_intersect_stages.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.cluster.err"
    benchmark:
        config["out_dir"] / "log/busco_intersect_stages.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.{busco_version}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        cpus=parameters["threads"]["busco5_intersect_stages"],
        node_options=parse_node_list("busco5_intersect_stages"),
        time=parameters["time"]["busco5_intersect_stages"],
        mem=parameters["memory_mb"]["busco5_intersect_stages"],
    threads:
        parameters["threads"]["busco5_intersect_stages"]
    shell:
         " OUTPUT_PREFIX={output.busco_legend};"
         " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.legend}};  "
         " workflow/scripts/busco/intersect_busco_results.py -b `echo '{input.busco_tables}' | tr ' ' ',' ` "
         "     -l {params.stages} -o ${{OUTPUT_PREFIX}} > {log.std} 2>&1;"
         " draw_features.py  -i {output.busco_counts_bedgraph} -o ${{OUTPUT_PREFIX}} -t bedgraph  "
         "     -n {output.busco_len} -z {output.busco_orderlist} -g {output.busco_legend} --hide_track_label "
         "     --color_column_name value -l {wildcards.busco_lineage} --figure_header_height 2 --subplots_adjust_top 0.7 "
         "     --x_tick_type int_number > {log.draw} 2>&1; "
         " draw_features.py  -i {output.no_complete_busco_counts_bedgraph} -o ${{OUTPUT_PREFIX}}.no_complete "
         "     -t bedgraph  -n {output.no_complete_busco_len} -z {output.busco_orderlist} -g {output.busco_legend} "
         "     --hide_track_label  --color_column_name value -l {wildcards.busco_lineage} "
         "     --figure_header_height 2 --subplots_adjust_top 0.7 "
         "     --x_tick_type int_number > {log.draw_no_complete} 2>&1; "
         " draw_features.py  -i {output.informative_busco_counts_bedgraph} -o ${{OUTPUT_PREFIX}}.informative "
         "     -t bedgraph  -n {output.informative_busco_len} -z {output.busco_orderlist} -g {output.busco_legend} "
         "     --hide_track_label  --color_column_name value -l {wildcards.busco_lineage} "
         "     --figure_header_height 2 --subplots_adjust_top 0.7 "
         "     --x_tick_type int_number > {log.draw_informative} 2>&1; "
         " touch {output.busco_counts_png}; "
         " touch {output.no_complete_busco_counts_png}; "
         " touch {output.informative_busco_counts_png}"


def get_labels_for_all_assemblies_in_chain(wildcards):
    chain_stage_dict = get_parameters_for_all_stages_in_chain(wildcards.parameters)
    for skip_stage in config["skip_qc_fasta_stage_list"]:
        if skip_stage in chain_stage_dict:
            chain_stage_dict.pop(skip_stage)
    label_list = []
    for stage in chain_stage_dict:
        for haplotype in stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"]:
            label_list.append("{0}_{1}".format(stage, haplotype))
    return label_list

def get_busco_tables_for_all_assemblies_in_chain(wildcards):
    chain_stage_dict = get_parameters_for_all_stages_in_chain(wildcards.parameters)
    for skip_stage in config["skip_qc_fasta_stage_list"]:
        if skip_stage in chain_stage_dict:
            chain_stage_dict.pop(skip_stage)
    busco_table_list = []
    for stage in chain_stage_dict:
        for haplotype in stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"]:
            busco_table_list += expand(config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/{genome_prefix}.{assembly_stage}.{haplotype}.{busco_lineage}.{busco_version}.full_table.tsv",
                                       assembly_stage=[stage],
                                       parameters=[chain_stage_dict[stage]],
                                       haplotype=[haplotype],
                                       allow_missing=True,)
    return busco_table_list


rule busco_intersect_all: # Downloading of busco datasets is performed by a different rule to avoid conflict between different instances of busco5
    priority: 500
    input:
        busco_tables=get_busco_tables_for_all_assemblies_in_chain,
        #config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.full_table.tsv",
    params:
        assemblies=lambda wildcards: ",".join(get_labels_for_all_assemblies_in_chain(wildcards))
    output:
        busco_legend=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/all_intersection/{genome_prefix}.{busco_lineage}.{busco_version}.legend",
        busco_orderlist=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/all_intersection/{genome_prefix}.{busco_lineage}.{busco_version}.orderlist",
        busco_merged_tsv=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/all_intersection/{genome_prefix}.{busco_lineage}.{busco_version}.merged.tsv",
        busco_len=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/all_intersection/{genome_prefix}.{busco_lineage}.{busco_version}.len",
        busco_counts_bedgraph=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/all_intersection/{genome_prefix}.{busco_lineage}.{busco_version}.counts.bedgraph",
        busco_counts_png=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/all_intersection/{genome_prefix}.{busco_lineage}.{busco_version}.png",
        no_complete_busco_len=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/all_intersection/{genome_prefix}.{busco_lineage}.{busco_version}.no_complete.len",
        no_complete_busco_counts_bedgraph=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/all_intersection/{genome_prefix}.{busco_lineage}.{busco_version}.no_complete.counts.bedgraph",
        no_complete_busco_counts_png=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/all_intersection/{genome_prefix}.{busco_lineage}.{busco_version}.no_complete.png",
        informative_busco_len=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/all_intersection/{genome_prefix}.{busco_lineage}.{busco_version}.informative.len",
        informative_busco_counts_bedgraph=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/all_intersection/{genome_prefix}.{busco_lineage}.{busco_version}.informative.counts.bedgraph",
        informative_busco_counts_png=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/{busco_version}/all_intersection/{genome_prefix}.{busco_lineage}.{busco_version}.informative.png",
    log:
        std=config["out_dir"] / "log/busco_intersect_all.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.{busco_version}.log",
        draw=config["out_dir"] / "log/busco_intersect_all.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.{busco_version}.draw.log",
        draw_no_complete=config["out_dir"] / "log/busco_intersect_all.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.{busco_version}.draw_no_complete.log",
        draw_informative=config["out_dir"] / "log/busco_intersect_all.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.{busco_version}.draw_informative.log",
        cluster_log=config["out_dir"] / "log/busco_intersect_all.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.{busco_version}.cluster.log",
        cluster_err=config["out_dir"] / "log/busco_intersect_all.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.{busco_version}.cluster.err"
    benchmark:
        config["out_dir"] / "log/busco_intersect_all.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.{busco_version}.benchmark.txt"
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
         " OUTPUT_PREFIX={output.busco_legend};"
         " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.legend}};  "
         " workflow/scripts/busco/intersect_busco_results.py -b `echo '{input.busco_tables}' | tr ' ' ',' ` "
         " -l {params.assemblies} -o ${{OUTPUT_PREFIX}} > {log.std} 2>&1;"
         " draw_features.py  -i {output.busco_counts_bedgraph} -o ${{OUTPUT_PREFIX}} -t bedgraph  "
         " -n {output.busco_len} -z {output.busco_orderlist} -g {output.busco_legend} --hide_track_label "
         " --color_column_name value -l {wildcards.busco_lineage} --figure_header_height 2 --subplots_adjust_top 0.7 "
         " --x_tick_type int_number > {log.draw} 2>&1; "
         " draw_features.py  -i {output.no_complete_busco_counts_bedgraph} -o ${{OUTPUT_PREFIX}}.no_complete "
         " -t bedgraph  -n {output.no_complete_busco_len} -z {output.busco_orderlist} -g {output.busco_legend} "
         " --hide_track_label  --color_column_name value -l {wildcards.busco_lineage} "
         " --figure_header_height 2 --subplots_adjust_top 0.7 "
         " --x_tick_type int_number > {log.draw_no_complete} 2>&1; "
         " draw_features.py  -i {output.informative_busco_counts_bedgraph} -o ${{OUTPUT_PREFIX}}.informative "
         " -t bedgraph  -n {output.informative_busco_len} -z {output.busco_orderlist} -g {output.busco_legend} "
         " --hide_track_label  --color_column_name value -l {wildcards.busco_lineage} "
         " --figure_header_height 2 --subplots_adjust_top 0.7 "
         " --x_tick_type int_number > {log.draw_informative} 2>&1; "
         " touch {output.busco_counts_png}; "
         " touch {output.no_complete_busco_counts_png}; "
         " touch {output.informative_busco_counts_png}"


rule create_busco_tracks:
    priority: 500
    input:
        busco_table="{fasta_dir}/assembly_qc/{busco_version}/{fasta_prefix}.{busco_lineage}.{busco_version}.full_table.tsv",
        log_dir=ancient("{fasta_dir}/log/")
    output:
        single_copy_track="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.{busco_lineage}.{busco_version}.single_copy.track.bed",
        duplicated_track="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.{busco_lineage}.{busco_version}.duplicated.track.bed",
        fragmented_track="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.{busco_lineage}.{busco_version}.fragmented.track.bed",
        single_copy_bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.{busco_lineage}.{busco_version}.single_copy.track.bedgraph",
        duplicated_bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.{busco_lineage}.{busco_version}.duplicated.track.bedgraph",
        fragmented_bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.{busco_lineage}.{busco_version}.fragmented.track.bedgraph",
    log:
        single_copy_track="{fasta_dir}/log/create_busco_tracks.{fasta_prefix}.{busco_lineage}.{busco_version}.single_copy.log",
        single_copy_grep="{fasta_dir}/log/create_busco_tracks.{fasta_prefix}.{busco_lineage}.{busco_version}.single_copy.grep.log",
        single_copy_awk="{fasta_dir}/log/create_busco_tracks.{fasta_prefix}.{busco_lineage}.{busco_version}.single_copy.awk.log",
        duplicated_track="{fasta_dir}/log/create_busco_tracks.{fasta_prefix}.{busco_lineage}.{busco_version}.duplicated.log",
        duplicated_grep="{fasta_dir}/log/create_busco_tracks.{fasta_prefix}.{busco_lineage}.{busco_version}.duplicated.grep.log",
        duplicated_awk="{fasta_dir}/log/create_busco_tracks.{fasta_prefix}.{busco_lineage}.{busco_version}.duplicated.awk.log",
        fragmented_track="{fasta_dir}/log/create_busco_tracks.{fasta_prefix}.{busco_lineage}.{busco_version}.fragmented.log",
        fragmented_grep="{fasta_dir}/log/create_busco_tracks.{fasta_prefix}.{busco_lineage}.{busco_version}.fragmented.grep.log",
        fragmented_awk="{fasta_dir}/log/create_busco_tracks.{fasta_prefix}.{busco_lineage}.{busco_version}.fragmented.awk.log",
        cluster_log="{fasta_dir}/log/create_busco_tracks.{fasta_prefix}.{busco_lineage}.{busco_version}.cluster.log",
        cluster_err="{fasta_dir}/log/create_busco_tracks.{fasta_prefix}.{busco_lineage}.{busco_version}.cluster.err"
    benchmark:
        "{fasta_dir}/log/create_busco_tracks.{fasta_prefix}.{busco_lineage}.{busco_version}.benchmark.txt"
    conda:
        config["conda"]["busco5.8"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco5.8"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("busco5_intersect_all"),
        cpus=parameters["threads"]["busco5_intersect_all"],
        time=parameters["time"]["busco5_intersect_all"],
        mem=parameters["memory_mb"]["busco5_intersect_all"],
    threads:
        parameters["threads"]["busco5_intersect_all"]
    shell:
        " grep -P '\\tComplete\\t' {input.busco_table} 2>{log.single_copy_grep} | awk -F '\\t' '{{ if ($6 == \"+\") {{printf \"%s\\t%i\\t%i\\t%s\\t%s\\n\",$3,$4,$5,$6,$1}} else {{printf \"%s\\t%i\\t%i\\t%s\\t%s\\n\",$3,$5,$4,$6,$1}} }}' 2>{log.single_copy_awk} | sort -k1,1V -k2,2n -k3,3n  > {output.single_copy_track} 2>{log.single_copy_track}; "
        " grep -P '\\tDuplicated\\t' {input.busco_table} 2>{log.duplicated_grep} | awk -F '\\t' '{{ if ($6 == \"+\") {{printf \"%s\\t%i\\t%i\\t%s\\t%s\\n\",$3,$4,$5,$6,$1}} else {{printf \"%s\\t%i\\t%i\\t%s\\t%s\\n\",$3,$5,$4,$6,$1}} }}' 2>{log.duplicated_awk} | sort -k1,1V -k2,2n -k3,3n  > {output.duplicated_track} 2>{log.duplicated_track}; "
        " grep -P '\\tFragmented\\t' {input.busco_table} 2>{log.fragmented_grep} | awk -F '\\t' '{{ if ($6 == \"+\") {{printf \"%s\\t%i\\t%i\\t%s\\t%s\\n\",$3,$4,$5,$6,$1}} else {{printf \"%s\\t%i\\t%i\\t%s\\t%s\\n\",$3,$5,$4,$6,$1}} }}' 2>{log.fragmented_awk} | sort -k1,1V -k2,2n -k3,3n  > {output.fragmented_track} 2>{log.fragmented_track}; "
        " awk -F'\\t' '{{printf \"%s\\t%i\\t%i\\t1\\n\",$1,$2,$3}}' {output.single_copy_track} > {output.single_copy_bedgraph} 2>>{log.single_copy_track}; "
        " awk -F'\\t' '{{printf \"%s\\t%i\\t%i\\t1\\n\",$1,$2,$3}}' {output.duplicated_track} > {output.duplicated_bedgraph} 2>>{log.duplicated_track}; "
        " awk -F'\\t' '{{printf \"%s\\t%i\\t%i\\t1\\n\",$1,$2,$3}}' {output.fragmented_track} > {output.fragmented_bedgraph} 2>>{log.fragmented_track}; "

rule create_busco_tracks_for_combined_haplotype:
    priority: 500
    input:
        single_haplotype_tracks=lambda wildcards: expand(config["out_dir"] / ("%s/%s/assembly_qc/tracks/%s.%s.{haplotype}/%s.%s.{haplotype}.%s.%s.{busco_type}.track.bed" % (wildcards.assembly_stage,
                                                                                                                                                                                  wildcards.parameters,
                                                                                                                                                                                  wildcards.genome_prefix,
                                                                                                                                                                                  wildcards.assembly_stage,
                                                                                                                                                                                  wildcards.genome_prefix,
                                                                                                                                                                                  wildcards.assembly_stage,
                                                                                                                                                                                  wildcards.busco_lineage,
                                                                                                                                                                                  wildcards.busco_version)),



                                                           busco_type=["single_copy", "duplicated", "fragmented"],
                                                           haplotype=stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"],
                                                           allow_missing=True)
    output:
        single_copy_track=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{merged_haplotype}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.{busco_lineage}.{busco_version}.single_copy.track.bed",
        duplicated_track=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{merged_haplotype}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.{busco_lineage}.{busco_version}.duplicated.track.bed",
        fragmented_track=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{merged_haplotype}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.{busco_lineage}.{busco_version}.fragmented.track.bed",
        single_copy_bedgraph=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{merged_haplotype}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.{busco_lineage}.{busco_version}.single_copy.track.bedgraph",
        duplicated_bedgraph=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{merged_haplotype}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.{busco_lineage}.{busco_version}.duplicated.track.bedgraph",
        fragmented_bedgraph=config["out_dir"] / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{merged_haplotype}/{genome_prefix}.{assembly_stage}.{merged_haplotype}.{busco_lineage}.{busco_version}.fragmented.track.bedgraph",
    params:
        haplotype_list=lambda wildcards: stage_dict[wildcards.assembly_stage].parameters[wildcards.parameters]["haplotype_list"],
        dir_prefix=lambda wildcards: config["out_dir"] / ("%s/%s/assembly_qc/tracks/%s.%s" % (wildcards.assembly_stage,
                                                                                        wildcards.parameters,
                                                                                        wildcards.genome_prefix,
                                                                                        wildcards.assembly_stage)),
        track_prefix=lambda wildcards: f"{wildcards.genome_prefix}.{wildcards.assembly_stage}",
        suffix=lambda wildcards: f"{wildcards.busco_lineage}.busco5",
    log:
        log=config["out_dir"] / "log/create_busco_tracks_for_combined_haplotype.{assembly_stage}.{parameters}.{genome_prefix}.{merged_haplotype}.{busco_lineage}.{busco_version}.log",
        cluster_log=config["out_dir"] / "log/create_busco_tracks_for_combined_haplotype.{assembly_stage}.{parameters}.{genome_prefix}.{merged_haplotype}.{busco_lineage}.{busco_version}.cluster.log",
        cluster_err=config["out_dir"] / "log/create_busco_tracks_for_combined_haplotype.{assembly_stage}.{parameters}.{genome_prefix}.{merged_haplotype}.{busco_lineage}.{busco_version}.cluster.err"
    benchmark:
        config["out_dir"] / "log/busco5.{assembly_stage}.{parameters}.{genome_prefix}.{merged_haplotype}.{busco_lineage}.{busco_version}.benchmark.txt"
    conda:
        config["conda"]["busco5.8"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco5.8"]["yaml"])
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
        " for BUSCO_TYPE in single_copy duplicated fragmented; "
        "   do "
        "   echo \"Processing ${{BUSCO_TYPE}}...\" >> {log.log} 2>&1; "
        "   OUT_TRACK={params.dir_prefix}.{wildcards.merged_haplotype}/{params.track_prefix}.{wildcards.merged_haplotype}.{params.suffix}.${{BUSCO_TYPE}}.track.bed; "
        "   > ${{OUT_TRACK}}; "
        "   for HAP in {params.haplotype_list}; "
        "       do "
        "       echo \"Processing ${{HAP}}...\" >> {log.log} 2>&1; "
        "       HAP_TRACK={params.dir_prefix}.${{HAP}}/{params.track_prefix}.${{HAP}}.{params.suffix}.${{BUSCO_TYPE}}.track.bed; "
        "       sed 's/^/'${{HAP}}'\./' ${{HAP_TRACK}} >> ${{OUT_TRACK}} 2>>{log.log}; "
        "       done;"
        "   OUT_BEDGRAPH={params.dir_prefix}.{wildcards.merged_haplotype}/{params.track_prefix}.{wildcards.merged_haplotype}.{params.suffix}.${{BUSCO_TYPE}}.track.bedgraph; "
        "   echo \"Converting to bedgraph...\" >> {log.log} 2>&1; "
        "   awk -F'\\t' '{{printf \"%s\\t%i\\t%i\\t1\\n\",$1,$2,$3}}' ${{OUT_TRACK}} > ${{OUT_BEDGRAPH}} 2>>{log.log}; "
        "   done; "


