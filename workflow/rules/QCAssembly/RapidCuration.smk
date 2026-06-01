import pandas as pd

localrules: select_scaffolds # create_curation_input_files_for_scaffolds, create_curation_input_files_for_contigs,

ruleorder: create_links_for_reference > select_scaffolds

#ruleorder: create_curation_input_files_for_scaffolds  > ref_faidx
#ruleorder: create_curation_input_files_for_contigs  > ref_faidx
#ruleorder: create_curation_input_files_for_scaffolds > get_seq_len
#ruleorder: create_curation_input_files_for_contigs > get_seq_len

def get_hic_bed_file(wildcards):
    #print(stage_dict["curation"]["prev_stage"]
    #print(stage_dict["curation"]["prev_stage"]["parameters"][wildcards.prev_stage_parameters])
    phasing_kmer_length = stage_dict[stage_dict["curation"]["prev_stage"]]["parameters"][wildcards.prev_stage_parameters]["option_set"]["phasing_kmer_length"] if "phasing_kmer_length" in stage_dict[stage_dict["curation"]["prev_stage"]]["parameters"][wildcards.prev_stage_parameters]["option_set"] else config["phasing_kmer_length"]
    return out_dir_path / "{0}/{1}/{2}/alignment/{3}/{4}.{0}.{3}.{2}.rmdup.bed".format(stage_dict["curation"]["prev_stage"],
                                                                                       wildcards.prev_stage_parameters,
                                                                                       wildcards.haplotype,
                                                                                       phasing_kmer_length,
                                                                                       wildcards.genome_prefix)

"""
rule create_curation_input_files_for_contigs: #
    input:
        fasta=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.new_contigs.fasta" % (stage_dict["curation"]["prev_stage"],
                                                                                                               stage_dict["curation"]["prev_stage"])),
        len=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.new_contigs.len" % (stage_dict["curation"]["prev_stage"],
                                                                                                           stage_dict["curation"]["prev_stage"])),
        fai=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.new_contigs.fasta.fai" % (stage_dict["curation"]["prev_stage"],
                                                                                                                 stage_dict["curation"]["prev_stage"])),
        transfer_agp=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.transfer.agp" % (stage_dict["curation"]["prev_stage"],
                                                                                                                 stage_dict["curation"]["prev_stage"])),
        #bed=get_hic_bed_file if not config["skip_higlass"] else []
    output:
        fasta=out_dir_path / "curation/{prev_stage_parameters, [^/]+}..{curation_parameters, [^/]+}/{haplotype, [^.]+}/contigs/{genome_prefix, [^/]+}.input.{haplotype}.fasta",
        len=out_dir_path / "curation/{prev_stage_parameters, [^/]+}..{curation_parameters, [^/]+}/{haplotype, [^.]+}/contigs/{genome_prefix, [^/]+}.input.{haplotype}.len",
        fai=out_dir_path / "curation/{prev_stage_parameters, [^/]+}..{curation_parameters, [^/]+}/{haplotype, [^.]+}/contigs/{genome_prefix, [^/]+}.input.{haplotype}.fasta.fai",
        transfer_agp=out_dir_path / "curation/{prev_stage_parameters, [^/]+}..{curation_parameters, [^/]+}/{haplotype, [^.]+}/contigs/{genome_prefix, [^/]+}.input.{haplotype}.transfer.agp",
        #bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.hic.bed" if not config["skip_higlass"] else [],
    log:
        cp=output_dict["log"]  / "create_curation_input_files.{prev_stage_parameters}..{curation_parameters}.contigs.{genome_prefix}.{haplotype}.cp.log",
        cluster_log=output_dict["cluster_log"] / "create_curation_input_files.{prev_stage_parameters}..{curation_parameters}.contigs.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_curation_input_files.{prev_stage_parameters}..{curation_parameters}.contigs.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_curation_input_files.{prev_stage_parameters}..{curation_parameters}.contigs.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_curation_input_files_for_contigs"),
        cpus=parameters["threads"]["create_curation_input_files"],
        time=parameters["time"]["create_curation_input_files"],
        mem=parameters["memory_mb"]["create_curation_input_files"]
    threads: parameters["threads"]["create_curation_input_files"]

    shell:
        " sed 's/\s.*//' `realpath -s {input.fasta}` > {output.fasta} 2>{log.cp}; "
        " samtools faidx -o {output.fai} {output.fasta} >> {log.cp} 2>&1; "
        " cp -f `realpath -s {input.len}` {output.len} >> {log.cp} 2>&1; "
        " cp -f `realpath -s {input.transfer_agp}` {output.transfer_agp} >> {log.cp} 2>&1; "

rule create_curation_bed_input_file: # Added as separated rule to allow turning on and off higlass track. DO NOT MERGE this rule with create_curation_input_files
    input:
        bed=get_hic_bed_file
    output:
        bed=out_dir_path / "curation/{prev_stage_parameters, [^/]+}..{curation_parameters, [^/]+}/{haplotype, [^.]+}/{seq_type, [^/]+}/{genome_prefix, [^/]+}.input.{haplotype}.hic.bed"
    log:
        cp=output_dict["log"]  / "create_bed_input_file.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cp.log",
        cluster_log=output_dict["cluster_log"] / "create_bed_input_file.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_bed_input_file.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_bed_input_file.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_curation_bed_input_file"),
        cpus=parameters["threads"]["create_curation_input_files"],
        time=parameters["time"]["create_curation_input_files"],
        mem=parameters["memory_mb"]["create_curation_input_files"]
    threads: parameters["threads"]["create_curation_input_files"]

    shell:
        " cp -f `realpath -s {input.bed}` {output.bed} > {log.cp} 2>&1; "
"""
rule select_scaffolds: #
    input:
        len="{len_dir}/{len_prefix}.len"
    output:
        whitelist="{len_dir}/{len_prefix, [^/]+}.{scaffold_length, [^/]+}.whitelist",
        orderlist="{len_dir}/{len_prefix, [^/]+}.{scaffold_length, [^/]+}.orderlist",
    params:
        max_scaffolds=lambda wildcards: parameters["tool_options"]["select_scaffolds"][wildcards.scaffold_length]["max_scaffolds"],
        max_len=lambda wildcards: parameters["tool_options"]["select_scaffolds"][wildcards.scaffold_length]["max_len"],
        min_len=lambda wildcards: parameters["tool_options"]["select_scaffolds"][wildcards.scaffold_length]["min_len"]
    log:
        std="{len_dir}/select_scaffolds.{len_prefix}.{scaffold_length}.std.log",
        cluster_log="{len_dir}/select_scaffolds.{len_prefix}.{scaffold_length}.{scaffold_length}.cluster.log",
        cluster_err="{len_dir}/select_scaffolds.{len_prefix}.{scaffold_length}.{scaffold_length}.cluster.err"
    benchmark:
        "{len_dir}/select_scaffolds.{len_prefix}.{scaffold_length}.{scaffold_length}.benchmark.txt"
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("select_scaffolds"),
        cpus=parameters["threads"]["select_scaffolds"],
        time=parameters["time"]["select_scaffolds"],
        mem=parameters["memory_mb"]["select_scaffolds"]
    threads: parameters["threads"]["select_scaffolds"]
    run:
        length_df = pd.read_csv(input.len, sep='\t', header=None, index_col=0, names=["scaffold", "length"])
        whitelist_sr = pd.Series(length_df[length_df["length"] >= params.min_len].index)
        if params.max_len:
            whitelist_sr = pd.Series(length_df[(length_df["length"] >= params.min_len) & (length_df["length"] <= params.max_len) ].index)
        if len(whitelist_sr) > params.max_scaffolds:
            whitelist_sr = whitelist_sr.iloc[:params.max_scaffolds]
        whitelist_sr.to_csv(output.whitelist, header=False, index=False)
        whitelist_sr.to_csv(output.orderlist, header=False, index=False)

rule create_windows: #
    input:
        len=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
    output:
        bed=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/{genome_prefix}.{assembly_stage}.{haplotype}.win{window, [0-9]+}.step{step, [0-9]+}.windows.bed",
    log:
        makewin=output_dict["log"]  / "create_windows.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.win{window}.step{step}.makewin.log",
        cluster_log=output_dict["cluster_log"] / "create_windows.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.win{window}.step{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_windows.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.win{window}.step{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_windows.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.win{window}.step{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_windows"),
        cpus=parameters["threads"]["create_windows"],
        time=parameters["time"]["create_windows"],
        mem=parameters["memory_mb"]["create_windows"],
        create_windows=1,
    threads: parameters["threads"]["create_windows"]

    shell:
        " bedtools makewindows -g {input.len} -w {wildcards.window} -s {wildcards.step} > {output.bed} 2>{log.makewin}; "

rule create_bedgraph_track: #
    input:
        track_bed=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/{track_type}/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.track.bed",
        windows_bed=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.win{window}.step{step}.windows.bed"
    output:
        bedgraph=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/tracks/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^./]+}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.win{window, [0-9]+}.step{step, [0-9]+}.track.bedgraph"
    log:
        intersect=output_dict["log"]  / "create_bedgraph_track.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.intersect.log",
        awk=output_dict["log"]  / "create_bedgraph_track.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.awk.log",
        map=output_dict["log"]  / "create_bedgraph_track.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.map.log",
        cluster_log=output_dict["cluster_log"] / "create_bedgraph_track.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_bedgraph_track.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_bedgraph_track.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_bedgraph_track"),
        cpus=parameters["threads"]["create_bedgraph_track"],
        time=parameters["time"]["create_bedgraph_track"],
        mem=parameters["memory_mb"]["create_bedgraph_track"],
        create_windows=1,
    threads: parameters["threads"]["create_bedgraph_track"]

    shell:
        " bedtools intersect -wao -a {input.windows_bed}  -b {input.track_bed}  2>{log.intersect} | "
        " awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$NF}}' 2>{log.awk} | "
        " workflow/scripts/sum_bed.py -c 3 > {output.bedgraph} 2>{log.map} "

rule scale_create_bedgraph_track: #
    input:
        bedgraph=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.win{window}.step{step}.track.bedgraph",
        yahs_juicer_pre_log=out_dir_path / "hic_scaffolding/{parameters}/{haplotype}/scaffolding/{genome_prefix}.hic_scaffolding.{haplotype}.log"

    output:
        bedgraph=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/scaled_tracks/{haplotype, [^./]+}/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype}.{track_type}.win{window, [0-9]+}.step{step, [0-9]+}.scaled.track.bedgraph"
    log:
        std=output_dict["log"]  / "scale_create_bedgraph_track.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.std.log",
        grep=output_dict["log"]  / "scale_create_bedgraph_track.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.grep.log",
        sed=output_dict["log"]  / "scale_create_bedgraph_track.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.sed.log",
        cluster_log=output_dict["cluster_log"] / "scale_create_bedgraph_track.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "scale_create_bedgraph_track.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "scale_create_bedgraph_track.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("scale_create_bedgraph_track"),
        cpus=parameters["threads"]["scale_create_bedgraph_track"],
        time=parameters["time"]["scale_create_bedgraph_track"],
        mem=parameters["memory_mb"]["scale_create_bedgraph_track"],
    threads: parameters["threads"]["scale_create_bedgraph_track"]

    shell:
        " SCALE=`grep 'scale factor:' {input.yahs_juicer_pre_log} 2>{log.grep} | sed 's/.*scale factor: //' 2>{log.sed}`; "
        " awk -v scale=${{SCALE}} '{{ print $1,$2/scale,$3/scale,$4 }}' {input.bedgraph} > {output.bedgraph} 2>{log.std}; "
"""
rule liftover_contig_bedgraph: #
    input:
        bedgraph=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/contigs/{genome_prefix}.input.{haplotype}.{track_type}.win{window}.step{step}.track.bedgraph",
        transfer_agp=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/contigs/{genome_prefix}.input.{haplotype}.transfer.agp",
    output:
        bedgraph=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/contigs/{genome_prefix}.input.{haplotype}.{track_type, [^./]+}.win{window}.step{step}.track.assembly.bedgraph"
    params:
        sorting_mem=parameters["memory_mb"]["liftover_contig_bedgraph"]
    log:
        std=output_dict["log"]  / "create_bedgraph_track.{prev_stage_parameters}..{curation_parameters}.contigs.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.std.log",
        cluster_log=output_dict["cluster_log"] / "create_bedgraph_track.{prev_stage_parameters}..{curation_parameters}.contigs.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_bedgraph_track.{prev_stage_parameters}..{curation_parameters}.contigs.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_bedgraph_track.{prev_stage_parameters}..{curation_parameters}.contigs.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("liftover_contig_bedgraph"),
        cpus=parameters["threads"]["liftover_contig_bedgraph"],
        time=parameters["time"]["liftover_contig_bedgraph"],
        mem=parameters["memory_mb"]["liftover_contig_bedgraph"] + 500,
        create_windows=1,
    threads: parameters["threads"]["liftover_contig_bedgraph"]

    shell:
        " ./workflow/scripts/curation/convert_contig_bed_to_assembly_bed.py -c {input.bedgraph} -t {input.transfer_agp} | "
        " sort --parallel {threads} -S {params.sorting_mem}M -k1,1V -k2,2n -k3,3n > {output.bedgraph} 2>{log.std}; "
"""
#print(parameters["tool_options"]["assembly_qc"])
#def tratat(wildcards):
#    print(parameters["tool_options"]["assembly_qc"])
#    return parse_option_flag("normalize_by_len", parameters["tool_options"]["assembly_qc"][wildcards.track_type], "-n")

rule get_track_stats: #
    input:
        bedgraph=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.win{window}.step{step}.track.bedgraph"
    output:
        per_scaffold_stat=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/track_stats/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^.]+}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type, [^./]+}.win{window}.step{step}.track.per_scaffold.stat",
        all_stat=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/track_stats/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^.]+}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type, [^./]+}.win{window}.step{step}.track.stat",
        thresholds=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/track_stats/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^.]+}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type, [^./]+}.win{window}.step{step}.track.thresholds"
    params:
        normalization=lambda wildcards: parse_option_flag("normalize_by_len", parameters["tool_options"]["assembly_qc"][wildcards.track_type], "-n")
    log:
        std=output_dict["log"]  / "get_track_stats.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.std.log",
        cluster_log=output_dict["cluster_log"] / "get_track_stats.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "get_track_stats.{assembly_stage}.{parameters}.{track_type}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "get_track_stats.{assembly_stage}.{parameters}.{track_type}..{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("get_track_stats"),
        cpus=parameters["threads"]["get_track_stats"],
        time=parameters["time"]["get_track_stats"],
        mem=parameters["memory_mb"]["get_track_stats"],
    threads: parameters["threads"]["get_track_stats"]

    shell:
        " OUTPUT_PREFIX={output.all_stat}; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.stat}}; "
        " workflow/scripts/curation/get_track_stats.py -i {input.bedgraph} {params.normalization} "
        " -o ${{OUTPUT_PREFIX}}  > {log.std} 2>&1; "

rule draw_track: #
    input:
        whitelist=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.{scaffold_length}.whitelist",
        orderlist=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.{scaffold_length}.orderlist",
        len_file=out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
        bedgraph=out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/tracks/{genome_prefix}.{assembly_stage}.{haplotype}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type}.win{window}.step{step}.track.bedgraph",
        relative_thresholds=lambda wildcards: expand(rules.get_track_stats.output.thresholds,
                                                       assembly_stage=[wildcards.assembly_stage],
                                                       parameters=[wildcards.parameters],
                                                       haplotype=[wildcards.haplotype],
                                                       genome_prefix=[wildcards.genome_prefix],
                                                       track_type=[wildcards.track_type],
                                                       window=[wildcards.window],
                                                       step=[wildcards.step]) if wildcards.threshold_type == 'relative' else [],
        all_stat=lambda wildcards: expand(rules.get_track_stats.output.all_stat,
                                            assembly_stage=[wildcards.assembly_stage],
                                            parameters=[wildcards.parameters],
                                            haplotype=[wildcards.haplotype],
                                            genome_prefix=[wildcards.genome_prefix],
                                            track_type=[wildcards.track_type],
                                            window=[wildcards.window],
                                            step=[wildcards.step]) if wildcards.threshold_type == 'relative' else []
    output:
        png=out_dir_path / "{assembly_stage, [^/]+}/{parameters, [^/]+}/assembly_qc/trackplots/{genome_prefix, [^/]+}.{assembly_stage}.{haplotype, [^./]+}/{genome_prefix}.{assembly_stage}.{haplotype}.{track_type, [^/]+}.{scaffold_length, [^./]+}.win{window}.step{step}.{threshold_type, [^/]+}.png"
    params: # TODO: move parameters from "curation" to "qc" or something similar
        thresholds=lambda wildcards: parse_option("absolute_thresholds",
                                                  parameters["tool_options"]["assembly_qc"][wildcards.track_type],
                                                  #stage_dict["curation"]["parameters"][wildcards.parameters]["option_set"][wildcards.track_type],
                                                  "--density_thresholds",
                                                  expression=lambda s: ",".join(list(map(str, s))))
    log:
        draw=output_dict["log"]  / "draw_track.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{track_type}.{scaffold_length}.win{window}.step{step}.{threshold_type}.draw.log",
        cluster_log=output_dict["cluster_log"] / "draw_track.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{track_type}.{scaffold_length}.win{window}.step{step}.{threshold_type}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "draw_track.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{track_type}.{scaffold_length}.win{window}.step{step}.{threshold_type}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "draw_track.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{track_type}.{scaffold_length}.win{window}.step{step}.{threshold_type}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("draw_track"),
        cpus=parameters["threads"]["draw_track"],
        time=parameters["time"]["draw_track"],
        mem=parameters["memory_mb"]["draw_track"],
    threads: parameters["threads"]["draw_track"]

    shell:
        " PREFIX={output.png}; "
        " PREFIX=${{PREFIX%.png}}; "
        " if [ '{wildcards.threshold_type}' == 'absolute' ]; "
        " then "
        "    THRESHOLDS='{params.thresholds}'; "
        "    TITLE={wildcards.track_type}; "
        " else "
        "    THRESHOLDS=\" --density_thresholds `head -n 1 {input.relative_thresholds} | sed 's/\\n//'`\"; "
        "    MEDIAN=`tail -n +2 {input.all_stat} | cut -f 4`; "
        "    TITLE=\"{wildcards.track_type}(median ${{MEDIAN}})\"; "
        " fi; "
        " draw_variant_window_densities.py -i {input.bedgraph} -t bedgraph -o ${{PREFIX}} -l \"${{TITLE}}\" "
        " -w {wildcards.window} -s {wildcards.step} --density_multiplier 1 "
        " -a {input.whitelist} -n {input.len_file} -z {input.orderlist} ${{THRESHOLDS}} "
        " --hide_track_label --rounded --subplots_adjust_left 0.15 --feature_name {wildcards.track_type} > {log.draw} 2>&1; "
