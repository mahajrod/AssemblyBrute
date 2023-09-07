
ruleorder: minimap2_purge_dups_qc > minimap2_purge_dups_reads
ruleorder: get_purge_dups_read_stat_qc > get_purge_dups_read_stat

rule minimap2_purge_dups_qc:
    input:
        fastq=lambda wildcards: output_dict["data"] / "fastq/{0}/filtered/{1}{2}".format(stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"],
                                                                                         wildcards.fileprefix,
                                                                                         config["fastq_extension"]),
        reference=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{genome_prefix}.purge_dups.{haplotype}.fasta"
    output:
        paf=out_dir_path  / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype, [^.]+}/{genome_prefix}.{haplotype}.{fileprefix}.paf.gz"
    params:
        index_size=lambda wildcards: parse_option("index_size", parameters["tool_options"]["minimap2"][stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"]], " -I "),
        alignment_scheme=lambda wildcards: parse_option("alignment_scheme", parameters["tool_options"]["minimap2"][stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"]], " -x "),
    log:
        std=output_dict["log"]  / "minimap2_purge_dups_qc.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}.log",
        cluster_log=output_dict["cluster_log"] / "minimap2_purge_dups_qc.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}..cluster.log",
        cluster_err=output_dict["cluster_error"] / "minimap2_purge_dups_qc.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}..cluster.err"
    benchmark:
        output_dict["benchmark"]  / "minimap2_purge_dups_qc.{prev_stage_parameters}.{purge_dups_parameters}.{haplotype}.{genome_prefix}.{fileprefix}..benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("minimap2_purge_dups_qc"),
        cpus=parameters["threads"]["minimap2"] ,
        time=parameters["time"]["minimap2"],
        mem=parameters["memory_mb"]["minimap2"]
    threads: parameters["threads"]["minimap2"]

    shell:
        " minimap2 {params.alignment_scheme} {params.index_size} -t {threads}  {input.reference} "
        " {input.fastq} 2>{log.std} |  gzip -c - > {output.paf}; "

rule get_purge_dups_read_stat_qc:
    input:
        paf=lambda wildcards: expand(rules.minimap2_purge_dups_qc.output.paf,
                           genome_prefix=[config["genome_prefix"]],
                           fileprefix=input_file_prefix_dict[stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["option_set"]["datatype"]],
                           allow_missing=True),
        genomescope_report=output_dict["kmer"] / "{0}/filtered/genomescope/{1}.{0}.filtered.{2}.{3}.genomescope.parameters".format(config["final_kmer_datatype"],
                                                                                                                                   config["genome_prefix"],
                                                                                                                                   config["final_kmer_length"],
                                                                                                                                   config["final_kmer_counter"]),
        before_pbstat=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/PB.stat",
        before_pbbasecov=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/PB.base.cov",
        before_cutoffs=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/cutoffs"
    output:
        pbstat=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype, [^.]+}/PB.stat",
        pbbasecov=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype, [^.]+}/PB.base.cov",
        cutoffs=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype, [^.]+}/cutoffs",
        coverage_plot=out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype, [^.]+}/{haplotype}.before-after.comparison.coverage.png"
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
        output_dict["benchmark"]  / "minimap2_purge_dups_read_stat_qc.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("get_purge_dups_read_stat"),
        cpus=parameters["threads"]["get_purge_dups_read_stat"] ,
        time=parameters["time"]["get_purge_dups_read_stat"],
        mem=parameters["memory_mb"]["get_purge_dups_read_stat"]
    threads: parameters["threads"]["get_purge_dups_read_stat"]

    shell:
        " OUT_DIR=`dirname {output.pbbasecov}`;"
        " COV_PLOT={output.coverage_plot}; "
        " COV_UPPER_BOUNDARY=`awk 'NR==2 {{printf \"%.0f\", {params.cov_multiplicator} * $2}}' {input.genomescope_report}`;"
        " if [ '{params.calcuts_upper_threshold}' != 'None' ] ; then COV_UPPER_BOUNDARY={params.calcuts_upper_threshold}; fi; "
        " pbcstat -O ${{OUT_DIR}} {input.paf} 1>{log.pbstat} 2>&1; "
        " calcuts -d 1 {params.calcuts_lower_threshold} {params.calcuts_haploid_diploid_threshold}"
        " -u ${{COV_UPPER_BOUNDARY}} {output.pbstat} > {output.cutoffs} 2>{log.calcuts}; " #check parameters for calcuts
        " workflow/scripts/purge_dups/draw_purge_dups_plot_all_haplotypes.py -b {input.before_pbstat},{output.pbstat} "
        " -l before,after -c {input.before_cutoffs},{output.cutoffs} -e png,svg -o ${{COV_PLOT%.png}} > {log.png} 2>&1; "

rule get_purge_stat_haplotype_comparison:
    input:
        before_pbstat=lambda wildcards: expand(out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/PB.stat",
                                               haplotype=stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["haplotype_list"],
                                               allow_missing=True,),
        before_cutoffs=lambda wildcards: expand(out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/{haplotype}/cutoffs",
                             haplotype=stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["haplotype_list"],
                             allow_missing=True,),
        after_pbstat=lambda wildcards: expand(out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype}/PB.stat",
                             haplotype=stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["haplotype_list"],
                             allow_missing=True,),
        after_cutoffs=lambda wildcards: expand(out_dir_path /  "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/{haplotype}/cutoffs",
                             haplotype=stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["haplotype_list"],
                             allow_missing=True,),
    output:
        before_coverage_plot=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/before.comparison.coverage.png",
        after_coverage_plot=out_dir_path / "purge_dups/{prev_stage_parameters}..{purge_dups_parameters}/assembly_qc/purge_dups/after.comparison.coverage.png"
    params:
        #label_list=lambda wildcards: stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["haplotype_list"],
        label_string=lambda wildcards: ",".join(stage_dict["purge_dups"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.purge_dups_parameters]["haplotype_list"])
    log:
        before=output_dict["log"] / "get_purge_stat_haplotype_comparison.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.before.log",
        after=output_dict["log"] / "get_purge_stat_haplotype_comparison.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.after.log",
        calcuts=output_dict["log"]  / "get_purge_stat_haplotype_comparison.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.calcuts.log",
        cluster_log=output_dict["cluster_log"] / "get_purge_stat_haplotype_comparison.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.cluster.log",
        cluster_err=output_dict["cluster_error"] / "get_purge_stat_haplotype_comparison.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "get_purge_stat_haplotype_comparison.{prev_stage_parameters}.{purge_dups_parameters}.purge_dups.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("get_purge_stat_haplotype_comparison"),
        cpus=parameters["threads"]["get_purge_stat_haplotype_comparison"] ,
        time=parameters["time"]["get_purge_stat_haplotype_comparison"],
        mem=parameters["memory_mb"]["get_purge_stat_haplotype_comparison"]
    threads: parameters["threads"]["get_purge_stat_haplotype_comparison"]

    shell:
        " BEFORE_COV_PLOT={output.before_coverage_plot}; "
        " AFTER_COV_PLOT={output.after_coverage_plot}; "
        " workflow/scripts/purge_dups/draw_purge_dups_plot_all_haplotypes.py "
        " -b `echo {input.before_pbstat} | tr ' ' ','` "
        " -l {params.label_string} -c `echo {input.before_cutoffs} | tr ' ' ','` "
        " -e png,svg -o ${{BEFORE_COV_PLOT%.png}} > {log.before} 2>&1; "
        " workflow/scripts/purge_dups/draw_purge_dups_plot_all_haplotypes.py "
        " -b `echo {input.after_pbstat} | tr ' ' ','` "
        " -l {params.label_string} -c `echo {input.after_cutoffs} | tr ' ' ','` "
        " -e png,svg -o ${{AFTER_COV_PLOT%.png}} > {log.after} 2>&1; "
