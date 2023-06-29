
rule last_index: #
    input:
        masked_fasta=rules.maskfasta.output.masked_fasta
    output:
        bck=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.YASS.R11.soft.bck",
        prj=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.YASS.R11.soft.prj",
        ssp=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.YASS.R11.soft.ssp",
        tis=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.YASS.R11.soft.tis",
        des=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.YASS.R11.soft.des",
        sds=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.YASS.R11.soft.sds",
        suf=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.YASS.R11.soft.suf",
    log:
        index=output_dict["log"]  / "last_index.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.index.log",
        cluster_log=output_dict["cluster_log"] / "last_index.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "last_index.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "last_index.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["last_index"] ,
        time=parameters["time"]["last_index"],
        mem=parameters["memory_mb"]["last_index"]
    threads: parameters["threads"]["last_index"]
    shell:
        " OUTPUT_PREFIX={output.bck}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.bck}}; "
        " lastdb  -P {threads} -c -u YASS -R11 ${{OUTPUT_PREFIX}} {input.masked_fasta} > {log.index} 2>&1; "

rule last_alignment: #
    input:
        database=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/input/{genome_prefix}.input.{target_haplotype}.YASS.R11.soft.bck",
        fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{query_haplotype}/input/{genome_prefix}.input.{query_haplotype}.softmasked.fasta",
    output:
        maf=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/input/{genome_prefix}.input.wga.{query_haplotype}.to.{target_haplotype}.YASS.R11.soft.maf.gz",
        tab=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/input/{genome_prefix}.input.wga.{query_haplotype}.to.{target_haplotype}.YASS.R11.soft.tab.gz",
    params:
        per_thread_mem=parameters["memory_mb"]["last_alignment_per_thread"],
    log:
        lastall=output_dict["log"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.lastall.log",
        tee=output_dict["log"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.tee.log",
        convert=output_dict["log"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.convert.log",
        pigz=output_dict["log"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.pigz.log",
        cluster_log=output_dict["cluster_log"] / "last_alignment.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "last_alignment.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "last_alignment.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["last_alignment"],
        time=parameters["time"]["last_alignment"],
        mem=parameters["memory_mb"]["last_alignment_per_thread"] * parameters["threads"]["last_alignment"]
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

rule draw_alignment: #
    input:
        tab=rules.last_alignment.output.tab,#out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/input/{genome_prefix}.input.wga.{query_haplotype}.to.{target_haplotype}.YASS.R11.soft.tab.gz",
        target_whitelist=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/input/{genome_prefix}.input.{target_haplotype}.whitelist",
        target_orderlist=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/input/{genome_prefix}.input.{target_haplotype}.orderlist",
        query_whitelist=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{query_haplotype}/input/{genome_prefix}.input.{query_haplotype}.whitelist",
        query_orderlist=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{query_haplotype}/input/{genome_prefix}.input.{query_haplotype}.orderlist",
    output:
        png=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/input/{genome_prefix}.input.wga.{query_haplotype}.to.{target_haplotype}.YASS.R11.soft.png",
        svg=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{target_haplotype}/input/{genome_prefix}.input.wga.{query_haplotype}.to.{target_haplotype}.YASS.R11.soft.svg",
    params:
        per_thread_mem=parameters["memory_mb"]["last_alignment_per_thread"],
    log:
        dotplot=output_dict["log"]  / "draw_alignment.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.dotplot.log",
        cluster_log=output_dict["cluster_log"] / "draw_alignment.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "draw_alignment.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "draw_alignment.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{query_haplotype}.to.{target_haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["draw_alignment"],
        time=parameters["time"]["draw_alignment"],
        mem=parameters["memory_mb"]["draw_alignment"]
    threads: parameters["threads"]["draw_alignment"]
    shell:
        " OUTPUT_PREFIX={output.png}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.png}}; "
        " dotplot_from_last_tab.py -i {input.tab} -w {input.target_whitelist} -x {input.query_whitelist} "
        " -u {input.target_orderlist} -z {input.query_orderlist} "
        " -l {wildcards.target_haplotype} -r {wildcards.query_haplotype} -e png,svg --axes_label_distance 7"
        " --bottom_offset 0.15 --top_offset 0.15 -o ${{OUTPUT_PREFIX}} > {log.dotplot} 2>&1; "