localrules: busco5_intersect_haplotypes, busco5_intersect_stages

rule busco5_download:
    priority: 500
    output:
        lineage_dir=directory(out_dir_path / "download/busco5/lineages/{busco_lineage}"),
    params:
        busco_download_dir=out_dir_path / "download/busco5/"
    log:
        std=output_dict["log"] / "busco5_download.{busco_lineage}.log",
        cluster_log=output_dict["cluster_log"] / "busco5_download.{busco_lineage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "busco5_download.{busco_lineage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "busco5_download.{busco_lineage}.benchmark.txt"
    conda:
        config["conda"]["busco"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco"]["yaml"])
    resources:
        cpus=parameters["threads"]["busco5_download"],
        time=parameters["time"]["busco5_download"],
        mem=parameters["memory_mb"]["busco5_download"],
    threads:
        parameters["threads"]["busco5_download"]
    shell:
         " busco --download_path {params.busco_download_dir} --download {wildcards.busco_lineage} > {log.std} 2>&1; "

rule busco5: # Downloading of busco datasets is performed by a different rule to avoid conflict between different instances of busco5
    priority: 500
    input:
        busco_lineage=rules.busco5_download.output.lineage_dir,
        assembly=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.fasta"
    output:
        tar_gz=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.tar.gz",
        summary=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.summary",
        summary_json=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.summary.json",
        busco_table=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.full_table.tsv",
        missing_busco_ids=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.missing.ids",
    log:
        std=output_dict["log"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.log",
        pigz=output_dict["log"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.pigz.log",
        cp=output_dict["log"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.cp.log",
        cluster_log=output_dict["cluster_log"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "busco5.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.benchmark.txt"
    conda:
        config["conda"]["busco"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco"]["yaml"])
    resources:
        cpus=parameters["threads"]["busco5"],
        time=parameters["time"]["busco5"],
        mem=parameters["memory_mb"]["busco5"],
    threads:
        parameters["threads"]["busco5"]
    shell:
         " BUSCO_DIR={output.tar_gz}; "
         " BUSCO_DIR=${{BUSCO_DIR%.tar.gz}}; "
         " busco --offline -m genome -l {input.busco_lineage} -c {threads} -i {input.assembly} "
         " -o `basename ${{BUSCO_DIR}}` --out_path `dirname ${{BUSCO_DIR}}` > {log.std} 2>&1;"
         " cp ${{BUSCO_DIR}}/short_summary.specific.{wildcards.busco_lineage}.{wildcards.genome_prefix}.{wildcards.assembly_stage}.{wildcards.haplotype}.busco5.{wildcards.busco_lineage}.txt {output.summary} ; "
         " cp ${{BUSCO_DIR}}/short_summary.specific.{wildcards.busco_lineage}.{wildcards.genome_prefix}.{wildcards.assembly_stage}.{wildcards.haplotype}.busco5.{wildcards.busco_lineage}.json {output.summary_json} ; "
         " cp ${{BUSCO_DIR}}/run_{wildcards.busco_lineage}/full_table.tsv {output.busco_table} ; "
         " cp ${{BUSCO_DIR}}/run_{wildcards.busco_lineage}/missing_busco_list.tsv {output.missing_busco_ids} ; "
         " tar cf - ${{BUSCO_DIR}} | pigz -p {threads} > {output.tar_gz} 2>{log.pigz} ;"
         " rm -r ${{BUSCO_DIR}}; "



rule busco5_intersect_haplotypes: # Downloading of busco datasets is performed by a different rule to avoid conflict between different instances of busco5
    priority: 500
    input:
        busco_tables=expand(rules.busco5.output.busco_table,
                           haplotype=haplotype_list,
                           allow_missing=True,)
        #out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.full_table.tsv",
    params:
        haplotypes=",".join(haplotype_list)
    output:
        busco_legend=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{busco_lineage}.busco.legend",
        busco_orderlist=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{busco_lineage}.busco.orderlist",
        busco_merged_tsv=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{busco_lineage}.busco.merged.tsv",
        busco_len=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{busco_lineage}.busco.len",
        busco_counts_bedgraph=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{busco_lineage}.busco.counts.bedgraph",
        no_complete_busco_len=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{busco_lineage}.no_complete.busco.len",
        no_complete_busco_counts_bedgraph=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{busco_lineage}.no_complete.busco.counts.bedgraph",
        informative_busco_len=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{busco_lineage}.informative.busco.len",
        informative_busco_counts_bedgraph=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{busco_lineage}.informative.busco.counts.bedgraph",
    log:
        std=output_dict["log"] / "busco5_intersect_haplotypes.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.log",
        pigz=output_dict["log"] / "busco5_intersect_haplotypes.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.pigz.log",
        cp=output_dict["log"] / "busco5_intersect_haplotypes.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.cp.log",
        cluster_log=output_dict["cluster_log"] / "busco5_intersect_haplotypes.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "busco5_intersect_haplotypes.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "busco5_intersect_haplotypes.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.benchmark.txt"
    conda:
        config["conda"]["busco"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco"]["yaml"])
    resources:
        cpus=parameters["threads"]["busco5_intersect_haplotypes"],
        time=parameters["time"]["busco5_intersect_haplotypes"],
        mem=parameters["memory_mb"]["busco5_intersect_haplotypes"],
    threads:
        parameters["threads"]["busco5_intersect_haplotypes"]
    shell:
         " OUTPUT_PREFIX={output.busco_legend};"
         " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.busco.legend}};  "
         " workflow/scripts/busco/intersect_busco_results.py -b `echo '{input.busco_tables}' | tr ' ' ',' ` "
         " -l {params.haplotypes} -o ${{OUTPUT_PREFIX}} > {log.std} 2>&1; "


rule busco5_intersect_stages:
    priority: 500
    input:
        busco_tables=lambda wildcards: [expand(rules.busco5.output.busco_table,
                                               assembly_stage=[stage],
                                               parameters=[get_parameters_for_all_stages_in_chain(wildcards.parameters)[stage]],
                                               #haplotype=haplotype_list,
                                               allow_missing=True,) for stage in get_parameters_for_all_stages_in_chain(wildcards.parameters)]
        #out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.full_table.tsv",
    params:
        stages=lambda wildcards: ",".join(get_parameters_for_all_stages_in_chain(wildcards.parameters))
    output:
        busco_legend=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.busco.legend",
        busco_orderlist=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.busco.orderlist",
        busco_merged_tsv=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.busco.merged.tsv",
        busco_len=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.busco.len",
        busco_counts_bedgraph=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.busco.counts.bedgraph",
        no_complete_busco_len=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.no_complete.busco.len",
        no_complete_busco_counts_bedgraph=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.no_complete.busco.counts.bedgraph",
        informative_busco_len=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.informative.busco.len",
        informative_busco_counts_bedgraph=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.informative.busco.counts.bedgraph",
    log:
        std=output_dict["log"] / "busco5_intersect_stages.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.log",
        pigz=output_dict["log"] / "busco5_intersect_stages.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.pigz.log",
        cp=output_dict["log"] / "busco5_intersect_stages.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.cp.log",
        cluster_log=output_dict["cluster_log"] / "busco5_intersect_stages.{assembly_stage}.{parameters}.{haplotype}.{genome_prefix}.{busco_lineage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "busco5_intersect_stages.{assembly_stage}.{parameters}.{haplotype}.{genome_prefix}.{busco_lineage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "busco5_intersect_stages.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{busco_lineage}.benchmark.txt"
    conda:
        config["conda"]["busco"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco"]["yaml"])
    resources:
        cpus=parameters["threads"]["busco5_intersect_stages"],
        time=parameters["time"]["busco5_intersect_stages"],
        mem=parameters["memory_mb"]["busco5_intersect_stages"],
    threads:
        parameters["threads"]["busco5_intersect_stages"]
    shell:
         " OUTPUT_PREFIX={output.busco_legend};"
         " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.busco.legend}};  "
         " workflow/scripts/busco/intersect_busco_results.py -b `echo '{input.busco_tables}' | tr ' ' ',' ` "
         " -l {params.stages} -o ${{OUTPUT_PREFIX}} > {log.std} 2>&1;"


def get_labels_for_all_assemblies_in_chain(wildcards):
    chain_stage_dict = get_parameters_for_all_stages_in_chain(wildcards.parameters)
    label_list = []
    for stage in chain_stage_dict:
        for haplotype in haplotype_list:
            label_list.append("{0}_{1}".format(stage, haplotype))
    return label_list

def get_busco_tables_for_all_assemblies_in_chain(wildcards):
    chain_stage_dict = get_parameters_for_all_stages_in_chain(wildcards.parameters)
    busco_table_list = []
    for stage in chain_stage_dict:
        for haplotype in haplotype_list:
            busco_table_list += expand(rules.busco5.output.busco_table,
                                       assembly_stage=[stage],
                                       haplotype=[haplotype],
                                       allow_missing=True,)
    return busco_table_list

rule busco5_intersect_all: # Downloading of busco datasets is performed by a different rule to avoid conflict between different instances of busco5
    priority: 500
    input:
        busco_tables=get_busco_tables_for_all_assemblies_in_chain,
        #out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype,[^.]+}.busco5.{busco_lineage}.full_table.tsv",
    params:
        assemblies=lambda wildcards: ",".join(get_labels_for_all_assemblies_in_chain(wildcards))
    output:
        busco_legend=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.busco.legend",
        busco_orderlist=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.busco.orderlist",
        busco_merged_tsv=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.busco.merged.tsv",
        busco_len=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.busco.len",
        busco_counts_bedgraph=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.busco.counts.bedgraph",
        no_complete_busco_len=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.no_complete.busco.len",
        no_complete_busco_counts_bedgraph=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.no_complete.busco.counts.bedgraph",
        informative_busco_len=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.informative.busco.len",
        informative_busco_counts_bedgraph=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.informative.busco.counts.bedgraph",
    log:
        std=output_dict["log"] / "busco5_intersect_all.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.log",
        pigz=output_dict["log"] / "busco5_intersect_all.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.pigz.log",
        cp=output_dict["log"] / "busco5_intersect_all.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.cp.log",
        cluster_log=output_dict["cluster_log"] / "busco5_intersect_all.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "busco5_intersect_all.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "busco5_intersect_all.{assembly_stage}.{parameters}.{genome_prefix}.{busco_lineage}.benchmark.txt"
    conda:
        config["conda"]["busco"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["busco"]["yaml"])
    resources:
        cpus=parameters["threads"]["busco5_intersect_all"],
        time=parameters["time"]["busco5_intersect_all"],
        mem=parameters["memory_mb"]["busco5_intersect_all"],
    threads:
        parameters["threads"]["busco5_intersect_all"]
    shell:
         " OUTPUT_PREFIX={output.busco_legend};"
         " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.busco.legend}};  "
         " workflow/scripts/busco/intersect_busco_results.py -b `echo '{input.busco_tables}' | tr ' ' ',' ` "
         " -l {params.assemblies} -o ${{OUTPUT_PREFIX}} > {log.std} 2>&1;"