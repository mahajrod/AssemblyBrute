"""
rule last_index: #
    input:
        masked_fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/scaffolds/{genome_prefix}.input.{haplotype}.softmasked.fasta"
    output:
        bck=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/scaffolds/{genome_prefix}.input.{haplotype}.YASS.R11.soft.bck",
        prj=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/scaffolds/{genome_prefix}.input.{haplotype}.YASS.R11.soft.prj",
        ssp=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/scaffolds/{genome_prefix}.input.{haplotype}.YASS.R11.soft.ssp",
        tis=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/scaffolds/{genome_prefix}.input.{haplotype}.YASS.R11.soft.tis",
        des=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/scaffolds/{genome_prefix}.input.{haplotype}.YASS.R11.soft.des",
        sds=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/scaffolds/{genome_prefix}.input.{haplotype}.YASS.R11.soft.sds",
        suf=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/scaffolds/{genome_prefix}.input.{haplotype}.YASS.R11.soft.suf",
    log:
        index=output_dict["log"]  / "last_index.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{haplotype}.index.log",
        cluster_log=output_dict["cluster_log"] / "last_index.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "last_index.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "last_index.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        cpus=parameters["threads"]["last_index"] ,
        time=parameters["time"]["last_index"],
        mem=parameters["memory_mb"]["last_index"]
    threads: parameters["threads"]["last_index"]
    shell:
        " OUTPUT_PREFIX={output.bck}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.bck}}; "
        " lastdb  -P {threads} -c -u YASS -R11 ${{OUTPUT_PREFIX}} {input.masked_fasta} > {log.index} 2>&1; "
"""

rule last_index: #
    input:
        masked_fasta="{fasta_prefix}.softmasked.fasta"
    output:
        bck="{fasta_prefix}.YASS.R11.soft.bck",
        prj="{fasta_prefix}.YASS.R11.soft.prj",
        ssp="{fasta_prefix}.YASS.R11.soft.ssp",
        tis="{fasta_prefix}.YASS.R11.soft.tis",
        des="{fasta_prefix}.YASS.R11.soft.des",
        sds="{fasta_prefix}.YASS.R11.soft.sds",
        suf="{fasta_prefix}.YASS.R11.soft.suf",
    log:
        index="{fasta_prefix}.index.log",
        cluster_log="{fasta_prefix}.cluster.log",
        cluster_err="{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("last_index"),
        cpus=parameters["threads"]["last_index"] ,
        time=parameters["time"]["last_index"],
        mem=parameters["memory_mb"]["last_index"]
    threads: parameters["threads"]["last_index"]
    shell:
        " OUTPUT_PREFIX={output.bck}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.bck}}; "
        " lastdb  -P {threads} -c -u YASS -R11 ${{OUTPUT_PREFIX}} {input.masked_fasta} > {log.index} 2>&1; "

def select_database(wildcards):
    if wildcards.target_haplotype in stage_dict["curation"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.curation_parameters]["haplotype_list"]:
        return out_dir_path / "curation/{0}..{1}/{2}/scaffolds/{3}.input.{2}.YASS.R11.soft.bck".format(wildcards.prev_stage_parameters,
                                                                                                       wildcards.curation_parameters,
                                                                                                       wildcards.target_haplotype,
                                                                                                       wildcards.genome_prefix)
    elif wildcards.target_haplotype in input_reference_filedict:
        return out_dir_path / "data/reference/{0}/{0}.YASS.R11.soft.bck".format(wildcards.target_haplotype),
    else:
        raise ValueError("ERROR!!! Unrecognized target haplotype: {}".format(wildcards.target_haplotype))

def select_query(wildcards):
    if wildcards.query_haplotype in stage_dict["curation"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.curation_parameters]["haplotype_list"]:
        return out_dir_path / "curation/{0}..{1}/{2}/scaffolds/{3}.input.{2}.softmasked.fasta".format(wildcards.prev_stage_parameters,
                                                                                                       wildcards.curation_parameters,
                                                                                                       wildcards.query_haplotype,
                                                                                                       wildcards.genome_prefix)
    elif wildcards.query_haplotype in input_reference_filedict:
        return out_dir_path / "data/reference/{0}/{0}.softmasked.fasta".format(wildcards.query_haplotype)
    else:
        raise ValueError("ERROR!!! Unrecognized query haplotype: {}".format(wildcards.query_haplotype))

rule last_alignment: #
    input:
        database=select_database,#out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/scaffolds/{genome_prefix}.input.{target_haplotype}.YASS.R11.soft.bck",
        fasta=select_query
    output:
        maf=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/scaffolds/{genome_prefix}.input.wga.{query_haplotype}.to.{target_haplotype}.YASS.R11.soft.min_len0.maf.gz",
        tab=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/scaffolds/{genome_prefix}.input.wga.{query_haplotype}.to.{target_haplotype}.YASS.R11.soft.min_len0.tab.gz",
    params:
        per_thread_mem=parameters["memory_mb"]["last_alignment_per_thread"],
    log:
        lastall=output_dict["log"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.lastall.log",
        tee=output_dict["log"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.tee.log",
        convert=output_dict["log"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.convert.log",
        pigz=output_dict["log"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.pigz.log",
        cluster_log=output_dict["cluster_log"] / "last_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "last_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("last_alignment"),
        cpus=parameters["threads"]["last_alignment"],
        time=parameters["time"]["last_alignment"],
        mem=parameters["memory_mb"]["last_alignment"]
    threads: parameters["threads"]["last_alignment"]
    shell:
        " INDEX_PREFIX={input.database}; "
        " INDEX_PREFIX=${{INDEX_PREFIX%.bck}}; "
        " MAF={output.maf}; "
        " MAF=${{MAF%.gz}}; "
        " TAB={output.tab}; "
        " TAB=${{TAB%.gz}}; "
        " lastal  -P {threads} -R11 -f MAF  ${{INDEX_PREFIX}} {input.fasta}  2>{log.lastall} | "
        " tee ${{MAF}} 2>{log.tee} | maf-convert tab > ${{TAB}} 2>{log.convert}; "
        " pigz -p {threads} ${{MAF}} ${{TAB}} > {log.pigz} 2>&1; "
"""
rule last_alignment_to_reference: #
    input:
        database=out_dir_path / "data/reference/{ref_name}/{ref_name}.YASS.R11.soft.bck",
        fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{query_haplotype}/scaffolds/{genome_prefix}.input.{query_haplotype}.softmasked.fasta",
    output:
        maf=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/scaffolds/{genome_prefix}.input.wga.{query_haplotype}.to.{ref_name}.YASS.R11.soft.min_len0.maf.gz",
        tab=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/scaffolds/{genome_prefix}.input.wga.{query_haplotype}.to.{ref_name}.YASS.R11.soft.min_len0.tab.gz",
    params:
        per_thread_mem=parameters["memory_mb"]["last_alignment_per_thread"],
    log:
        lastall=output_dict["log"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{ref_name}.lastall.log",
        tee=output_dict["log"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{ref_name}.tee.log",
        convert=output_dict["log"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{ref_name}.convert.log",
        pigz=output_dict["log"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{ref_name}.pigz.log",
        cluster_log=output_dict["cluster_log"] / "last_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{ref_name}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "last_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{ref_name}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{ref_name}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        cpus=parameters["threads"]["last_alignment"],
        time=parameters["time"]["last_alignment"],
        mem=parameters["memory_mb"]["last_alignment"]
    threads: parameters["threads"]["last_alignment"]
    shell:
        " INDEX_PREFIX={input.database}; "
        " INDEX_PREFIX=${{INDEX_PREFIX%.bck}}; "
        " MAF={output.maf}; "
        " MAF=${{MAF%.gz}}; "
        " TAB={output.tab}; "
        " TAB=${{TAB%.gz}}; "
        " lastal  -P {threads} -R11 -f MAF  ${{INDEX_PREFIX}} {input.fasta}  2>{log.lastall} | "
        " tee ${{MAF}} 2>{log.tee} | maf-convert tab > ${{TAB}} 2>{log.convert}; "
        " pigz -p {threads} ${{MAF}} ${{TAB}} > {log.pigz} 2>&1; "
"""
"""
rule filter_last_alignment_by_target_hit_len: #
    input:
        tab=rules.last_alignment.output.tab,
    output:
        tab=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/scaffolds/{genome_prefix}.input.wga.{query_haplotype}.to.{target_haplotype}.YASS.R11.soft.min_len{min_target_len}.tab.gz",
    params:
        per_thread_mem=parameters["memory_mb"]["last_alignment_per_thread"],
    log:
        zcat=output_dict["log"]  / "filter_last_alignment_by_len.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.min_len{min_target_len}.zcat.log",
        grep=output_dict["log"]  / "filter_last_alignment_by_len.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.min_len{min_target_len}.grep.log",
        awk=output_dict["log"]  / "filter_last_alignment_by_len.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.min_len{min_target_len}.awk.log",
        gzip=output_dict["log"]  / "filter_last_alignment_by_len.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.min_len{min_target_len}.gzip.log",
        cluster_log=output_dict["cluster_log"] / "filter_last_alignment_by_len.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.min_len{min_target_len}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "filter_last_alignment_by_len.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.min_len{min_target_len}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "filter_last_alignment_by_len.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.min_len{min_target_len}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        cpus=parameters["threads"]["filter_last_alignment_by_len"],
        time=parameters["time"]["filter_last_alignment_by_len"],
        mem=parameters["memory_mb"]["filter_last_alignment_by_len"]
    threads: parameters["threads"]["filter_last_alignment_by_len"]
    shell:
        " zcat {input.tab} 2>{log.zcat} | "
        " grep -vP '^#' 2>{log.grep} | "
        " awk -F'\t' '{{if ($4>={wildcards.min_target_len}) print $0}}' 2>{log.awk} | "
        " gzip -c > {output.tab} 2>{log.gzip}; "
"""
rule filter_last_alignment_by_target_hit_len: #
    input:
        tab="{tab_file_prefix}.YASS.R11.soft.min_len0.tab.gz"
    output:
        tab="{tab_file_prefix}.YASS.R11.soft.min_len{min_target_len}.tab.gz",
    params:
        per_thread_mem=parameters["memory_mb"]["last_alignment_per_thread"],
    log:
        zcat="{tab_file_prefix}.min_len{min_target_len}.zcat.log",
        grep="{tab_file_prefix}.min_len{min_target_len}.grep.log",
        awk="{tab_file_prefix}.min_len{min_target_len}.awk.log",
        gzip="{tab_file_prefix}.min_len{min_target_len}.gzip.log",
        cluster_log="{tab_file_prefix}.min_len{min_target_len}.cluster.log",
        cluster_err="{tab_file_prefix}.min_len{min_target_len}.cluster.err"
    benchmark:
        "{tab_file_prefix}.min_len{min_target_len}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("filter_last_alignment_by_target_hit_len"),
        cpus=parameters["threads"]["filter_last_alignment_by_len"],
        time=parameters["time"]["filter_last_alignment_by_len"],
        mem=parameters["memory_mb"]["filter_last_alignment_by_len"]
    threads: parameters["threads"]["filter_last_alignment_by_len"]
    shell:
        " zcat {input.tab} 2>{log.zcat} | "
        " grep -vP '^#' 2>{log.grep} | "
        " awk -F'\t' '{{if ($4>={wildcards.min_target_len}) print $0}}' 2>{log.awk} | "
        " gzip -c > {output.tab} 2>{log.gzip}; "

def select_query_whitelist(wildcards):
    if wildcards.query_haplotype in stage_dict["curation"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.curation_parameters]["haplotype_list"]:
        return out_dir_path / "curation/{0}..{1}/{2}/scaffolds/{3}.input.{2}.whitelist".format(wildcards.prev_stage_parameters,
                                                                                                       wildcards.curation_parameters,
                                                                                                       wildcards.query_haplotype,
                                                                                                       wildcards.genome_prefix)
    elif wildcards.query_haplotype in input_reference_filedict:
        return out_dir_path / "data/reference/{0}/{0}.whitelist".format(wildcards.query_haplotype)
    else:
        raise ValueError("ERROR!!! Unrecognized query haplotype: {}".format(wildcards.query_haplotype))

def select_query_orderlist(wildcards):
    if wildcards.query_haplotype in stage_dict["curation"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.curation_parameters]["haplotype_list"]:
        return out_dir_path / "curation/{0}..{1}/{2}/scaffolds/{3}.input.{2}.orderlist".format(wildcards.prev_stage_parameters,
                                                                                                       wildcards.curation_parameters,
                                                                                                       wildcards.query_haplotype,
                                                                                                       wildcards.genome_prefix)
    elif wildcards.query_haplotype in input_reference_filedict:
        return out_dir_path / "data/reference/{0}/{0}.orderlist".format(wildcards.query_haplotype)
    else:
        raise ValueError("ERROR!!! Unrecognized query haplotype: {}".format(wildcards.query_haplotype))

def select_query_synfile(wildcards):
    if wildcards.query_haplotype in stage_dict["curation"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.curation_parameters]["haplotype_list"]:
        return []
    elif wildcards.query_haplotype in input_reference_filedict:
        return out_dir_path / "data/reference/{0}/{0}.syn".format(wildcards.query_haplotype)
    else:
        raise ValueError("ERROR!!! Unrecognized query haplotype: {}".format(wildcards.query_haplotype))

rule draw_alignment: #
    input:
        tab=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/scaffolds/{genome_prefix}.input.wga.{query_haplotype}.to.{target_haplotype}.YASS.R11.soft.min_len{min_target_len}.tab.gz",#out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/input/{genome_prefix}.input.wga.{query_haplotype}.to.{target_haplotype}.YASS.R11.soft.tab.gz",
        target_whitelist=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/scaffolds/{genome_prefix}.input.{target_haplotype}.whitelist",
        target_orderlist=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/scaffolds/{genome_prefix}.input.{target_haplotype}.orderlist",
        query_whitelist=select_query_whitelist, #out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{query_haplotype}/scaffolds/{genome_prefix}.input.{query_haplotype}.whitelist",
        query_orderlist=select_query_orderlist,
        query_synfile=select_query_synfile,
    output:
        png=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/scaffolds/{genome_prefix}.input.wga.{query_haplotype}.to.{target_haplotype}.YASS.R11.soft.min_len{min_target_len}.png",
        svg=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/scaffolds/{genome_prefix}.input.wga.{query_haplotype}.to.{target_haplotype}.YASS.R11.soft.min_len{min_target_len}.svg",
    params:
        per_thread_mem=parameters["memory_mb"]["last_alignment_per_thread"],
        query_syn_file=lambda wildcards: " --query_syn_file {0}".format(str(select_query_synfile(wildcards))) if select_query_synfile(wildcards) else ""
    log:
        dotplot=output_dict["log"]  / "draw_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.min_len{min_target_len}.dotplot.log",
        cluster_log=output_dict["cluster_log"] / "draw_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.min_len{min_target_len}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "draw_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.min_len{min_target_len}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "draw_alignment.{prev_stage_parameters}..{curation_parameters}.scaffolds.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.min_len{min_target_len}.benchmark.txt"
    conda:
        config["conda"]["chromodoter"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["chromodoter"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("draw_alignment"),
        cpus=parameters["threads"]["draw_alignment"],
        time=parameters["time"]["draw_alignment"],
        mem=parameters["memory_mb"]["draw_alignment"]
    threads: parameters["threads"]["draw_alignment"]
    shell:
        " OUTPUT_PREFIX={output.png}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.png}}; "
        " dotplot_from_last_tab.py -i {input.tab} -w {input.target_whitelist} -x {input.query_whitelist} "
        " -u {input.target_orderlist} -z {input.query_orderlist} "
        " {params.query_syn_file} --query_syn_file_key_column 1 --query_syn_file_value_column 0"
        " -l {wildcards.target_haplotype} -r {wildcards.query_haplotype} -e png,svg --axes_label_distance 6 "
        " --bottom_offset 0.15 --left_offset 0.15 --top_offset 0.9 -o ${{OUTPUT_PREFIX}} > {log.dotplot} 2>&1; "