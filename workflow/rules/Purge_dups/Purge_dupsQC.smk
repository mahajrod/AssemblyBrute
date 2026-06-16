
ruleorder: minimap2_purge_dups_qc > minimap2_purge_dups_reads
ruleorder: get_purge_dups_read_stat_qc > get_purge_dups_read_stat


rule minimap2_purge_dups_qc:
    input:
        fastq=lambda wildcards: output_dict["data"] / "fastq/{0}/filtered/{1}{2}".format(wildcards.datatype, #stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"],
                                                                                         wildcards.fileprefix,
                                                                                         config["fastq_extension"]),
        reference=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{genome_prefix}.purge_dups.{haplotype}.fasta"
    output:
        paf=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype}/{datatype}/{genome_prefix}.{haplotype}.{fileprefix}.paf.gz"
    params:
        index_size=lambda wildcards: parse_option("index_size", parameters["tool_options"]["minimap2"][wildcards.datatype], " -I "),
        alignment_scheme=lambda wildcards: parse_option("alignment_scheme", parameters["tool_options"]["minimap2"][wildcards.datatype], " -x "),
    log:
        std=output_dict["log"]  / "minimap2_purge_dups_qc.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{datatype}.{fileprefix}.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_purge_dups_qc.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{datatype}.{genome_prefix}.{fileprefix}..cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_purge_dups_qc.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{datatype}.{genome_prefix}.{fileprefix}..cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_qc.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{datatype}.{fileprefix}..benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("minimap2_purge_dups_qc"),
        cpus=parameters["threads"]["minimap2"] ,
        time=parameters["time"]["minimap2"],
        mem=parameters["memory_mb"]["minimap2"]
    threads: parameters["threads"]["minimap2"]

    shell:
        " minimap2 {params.alignment_scheme} {params.index_size} -t {threads}  {input.reference} "
        " {input.fastq} 2>{log.std} |  gzip -c - > {output.paf}; "

def get_paf_list_for_qc(wildcards):
    paf_list = []
    for datatype in stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["main_datatypes"]:
        paf_list += expand(rules.minimap2_purge_dups_qc.output.paf,
                           datatype=[datatype],
                           fileprefix=input_file_prefix_dict[datatype],
                           genome_prefix=[config["genome_prefix"]],
                           allow_missing=True)
    return paf_list

rule get_purge_dups_read_stat_qc:
    input:
        paf=get_paf_list_for_qc,
        #paf=lambda wildcards: expand(rules.minimap2_purge_dups_qc.output.paf,
        #                   genome_prefix=[config["genome_prefix"]],
        #                   fileprefix=input_file_prefix_dict[stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"]],
        #                   allow_missing=True),
        genomescope_report=output_dict["kmer"] / "{0}/filtered/genomescope/{1}.{0}.filtered.{2}.{3}.genomescope.parameters".format("_".join(config["final_kmer_datatypes"]),
                                                                                                                                   config["genome_prefix"],
                                                                                                                                   config["final_kmer_length"],
                                                                                                                                   config["final_kmer_counter"]),

    output:
        pbstat=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype}/PB.stat",
        pbbasecov=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype}/PB.base.cov",
        cutoffs=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype}/cutoffs",
    params:
        cov_multiplicator=lambda wildcards: stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["cov_multiplicator"],
        calcuts_lower_threshold=lambda wildcards: parse_option("lower_threshold", config["tool_manually_adjusted_features"]["calcuts"], " -l "),
        calcuts_haploid_diploid_threshold=lambda wildcards: parse_option("haploid_diploid_threshold", config["tool_manually_adjusted_features"]["calcuts"], " -m "),
        calcuts_upper_threshold=str(config["tool_manually_adjusted_features"]["calcuts"]["upper_threshold"]), # None needs to be converted to "None"
    log:
        pbstat=output_dict["log"] / "get_purge_dups_read_stat_qc.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.pbstat.log",
        png=output_dict["log"] / "get_purge_dups_read_stat_qc.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.png.log",
        calcuts=output_dict["log"]  / "get_purge_dups_read_stat_qc.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.calcuts.log",
        cluster_log=output_dict["cluster_log"] / "get_purge_dups_read_stat_qc.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "get_purge_dups_read_stat_qc.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "get_purge_dups_read_stat_qc..{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("get_purge_dups_read_stat"),
        cpus=parameters["threads"]["get_purge_dups_read_stat"] ,
        time=parameters["time"]["get_purge_dups_read_stat"],
        mem=parameters["memory_mb"]["get_purge_dups_read_stat"]
    threads: parameters["threads"]["get_purge_dups_read_stat"]

    shell:
        " OUT_DIR=`dirname {output.pbbasecov}`;"
        " COV_UPPER_BOUNDARY=`awk 'NR==2 {{printf \"%.0f\", {params.cov_multiplicator} * $2}}' {input.genomescope_report}`;"
        " if [ '{params.calcuts_upper_threshold}' != 'None' ] ; then COV_UPPER_BOUNDARY={params.calcuts_upper_threshold}; fi; "
        " pbcstat -O ${{OUT_DIR}} {input.paf} 1>{log.pbstat} 2>&1; "
        " calcuts -d 1 {params.calcuts_lower_threshold} {params.calcuts_haploid_diploid_threshold}"
        " -u ${{COV_UPPER_BOUNDARY}} {output.pbstat} > {output.cutoffs} 2>{log.calcuts}; " #check parameters for calcuts


rule draw_before_after_plot:
    input:
        before_pbstat=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/first_stage/{haplotype}/PB.stat",
        before_pbbasecov=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/first_stage/{haplotype}/PB.base.cov",
        before_cutoffs=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/first_stage/{haplotype}/cutoffs",
        after_cutoffs=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype}/cutoffs",
        after_pbstat=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype}/PB.stat",
    output:
        coverage_plot=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype}/{haplotype}.before-after.comparison.coverage.png"

    log:
        pbstat=output_dict["log"] / "draw_before_after_plot.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.pbstat.log",
        png=output_dict["log"] / "draw_before_after_plot.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.png.log",
        calcuts=output_dict["log"]  / "gdraw_before_after_plot.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.calcuts.log",
        cluster_log=output_dict["cluster_log"] / "draw_before_after_plot.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "draw_before_after_plot.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "draw_before_after_plot{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("get_purge_dups_read_stat"),
        cpus=parameters["threads"]["get_purge_dups_read_stat"] ,
        time=parameters["time"]["get_purge_dups_read_stat"],
        mem=parameters["memory_mb"]["get_purge_dups_read_stat"]
    threads: parameters["threads"]["get_purge_dups_read_stat"]

    shell:
        " COV_PLOT={output.coverage_plot}; "
        " workflow/scripts/purge_dups/draw_purge_dups_plot_all_haplotypes.py -b {input.before_pbstat},{input.after_pbstat} "
        " -l before,after -c {input.before_cutoffs},{input.after_cutoffs} -e png,svg -o ${{COV_PLOT%.png}} > {log.png} 2>&1; "


rule draw_before_plot:
    input:
        before_pbstat=lambda wildcards: expand(out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/first_stage/{haplotype}/PB.stat",
                                               haplotype=stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["haplotype_list"],
                                               allow_missing=True,),
        before_cutoffs=lambda wildcards: expand(out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/first_stage/{haplotype}/cutoffs",
                             haplotype=stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["haplotype_list"],
                             allow_missing=True,),
    output:
        before_coverage_plot=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/before.comparison.coverage.png",
    params:
        #label_list=lambda wildcards: stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["haplotype_list"],
        label_string=lambda wildcards: ",".join(stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["haplotype_list"])
    log:
        before=output_dict["log"] / "draw_before_plot.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.before.log",
        cluster_log=output_dict["cluster_log"] / "draw_before_plot.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.cluster.log",
        cluster_err=output_dict["cluster_error"] / "draw_before_plot.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "draw_before_plot.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("get_purge_stat_haplotype_comparison"),
        cpus=parameters["threads"]["get_purge_stat_haplotype_comparison"] ,
        time=parameters["time"]["get_purge_stat_haplotype_comparison"],
        mem=parameters["memory_mb"]["get_purge_stat_haplotype_comparison"]
    threads: parameters["threads"]["get_purge_stat_haplotype_comparison"]

    shell:
        " BEFORE_COV_PLOT={output.before_coverage_plot}; "
        " workflow/scripts/purge_dups/draw_purge_dups_plot_all_haplotypes.py "
        " -b `echo {input.before_pbstat} | tr ' ' ','` "
        " -l {params.label_string} -c `echo {input.before_cutoffs} | tr ' ' ','` "
        " -e png,svg -o ${{BEFORE_COV_PLOT%.png}} > {log.before} 2>&1; "


rule draw_after_plot:
    input:
        after_pbstat=lambda wildcards: expand(out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype}/PB.stat",
                             haplotype=stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["haplotype_list"],
                             allow_missing=True,),
        after_cutoffs=lambda wildcards: expand(out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype}/cutoffs",
                             haplotype=stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["haplotype_list"],
                             allow_missing=True,),
    output:
        after_coverage_plot=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/after.comparison.coverage.png"
    params:
        label_string=lambda wildcards: ",".join(stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["haplotype_list"])
    log:
        after=output_dict["log"] / "draw_after_plot.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.after.log",
        cluster_log=output_dict["cluster_log"] / "draw_after_plot.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.cluster.log",
        cluster_err=output_dict["cluster_error"] / "draw_after_plot.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "draw_after_plot.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("get_purge_stat_haplotype_comparison"),
        cpus=parameters["threads"]["get_purge_stat_haplotype_comparison"] ,
        time=parameters["time"]["get_purge_stat_haplotype_comparison"],
        mem=parameters["memory_mb"]["get_purge_stat_haplotype_comparison"]
    threads: parameters["threads"]["get_purge_stat_haplotype_comparison"]

    shell:
        " AFTER_COV_PLOT={output.after_coverage_plot}; "
        " workflow/scripts/purge_dups/draw_purge_dups_plot_all_haplotypes.py "
        " -b `echo {input.after_pbstat} | tr ' ' ','` "
        " -l {params.label_string} -c `echo {input.after_cutoffs} | tr ' ' ','` "
        " -e png,svg -o ${{AFTER_COV_PLOT%.png}} > {log.after} 2>&1; "