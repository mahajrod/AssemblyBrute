import pandas as pd

localrules: create_curation_input_links, select_long_scaffolds, create_windows

def get_hic_bed_file(wildcards):
    #print(stage_dict["curation"]["prev_stage"]
    #print(stage_dict["curation"]["prev_stage"]["parameters"][wildcards.prev_stage_parameters])
    phasing_kmer_length = stage_dict[stage_dict["curation"]["prev_stage"]]["parameters"][wildcards.prev_stage_parameters]["option_set"]["phasing_kmer_length"] if "phasing_kmer_length" in stage_dict[stage_dict["curation"]["prev_stage"]]["parameters"][wildcards.prev_stage_parameters]["option_set"] else config["phasing_kmer_length"]
    return out_dir_path / "{0}/{1}/{2}/alignment/{3}/{4}.{0}.{3}.{2}.rmdup.bed".format(stage_dict["curation"]["prev_stage"],
                                                                                       wildcards.prev_stage_parameters,
                                                                                       wildcards.haplotype,
                                                                                       phasing_kmer_length,
                                                                                       wildcards.genome_prefix)

rule create_curation_input_links: #
    input:
        fasta=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta" % (stage_dict["curation"]["prev_stage"],
                                                                                                   stage_dict["curation"]["prev_stage"])),
        len=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.len" % (stage_dict["curation"]["prev_stage"],
                                                                                                   stage_dict["curation"]["prev_stage"])),
        fai=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta.fai" % (stage_dict["curation"]["prev_stage"],
                                                                                                       stage_dict["curation"]["prev_stage"])),
        bed=get_hic_bed_file
    output:
        fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.fasta",
        len=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.len",
        fai=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.fasta.fai",
        bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.hic.bed",
    log:
        ln=output_dict["log"]  / "create_curation_input_links.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "create_curation_input_links.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_curation_input_links.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_curation_input_links.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["create_curation_input_links"],
        time=parameters["time"]["create_curation_input_links"],
        mem=parameters["memory_mb"]["create_curation_input_links"]
    threads: parameters["threads"]["create_curation_input_links"]

    shell:
        " ln -sf `realpath -s {input.fasta}` {output.fasta} > {log.ln} 2>&1; "
        " ln -sf `realpath -s {input.fai}` {output.fai} >> {log.ln} 2>&1; "
        " ln -sf `realpath -s {input.bed}` {output.bed} >> {log.ln} 2>&1; "
        " ln -sf `realpath -s {input.len}` {output.len} >> {log.ln} 2>&1; "

rule select_long_scaffolds: #
    input:
        len=rules.create_curation_input_links.output.len
    output:
        whitelist=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.whitelist",
        orderlist=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.orderlist"
    log:
        ln=output_dict["log"]  / "select_long_scaffolds.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "select_long_scaffolds.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "select_long_scaffolds.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "select_long_scaffolds.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    #conda:
    #    config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["select_long_scaffolds"],
        time=parameters["time"]["select_long_scaffolds"],
        mem=parameters["memory_mb"]["select_long_scaffolds"]
    threads: parameters["threads"]["select_long_scaffolds"]
    run:
        length_df = pd.read_csv(input.len, sep='\t', header=None, index_col=0, names=["scaffold", "length"])
        threshold = length_df["length"].iloc[0] / 50
        whitelist_sr = pd.Series(length_df[length_df["length"] >= threshold].index)
        whitelist_sr.to_csv(output.whitelist, header=False, index=False)
        whitelist_sr.to_csv(output.orderlist, header=False, index=False)

rule create_windows: #
    input:
        len=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.len",
    output:
        bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.win{window}.step{step}.windows.bed",
    log:
        makewin=output_dict["log"]  / "create_windows.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.win{window}.step{step}.makewin.log",
        cluster_log=output_dict["cluster_log"] / "create_windows.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.win{window}.step{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_windows.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.win{window}.step{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_windows{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.win{window}.step{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["create_windows"],
        time=parameters["time"]["create_windows"],
        mem=parameters["memory_mb"]["create_windows"],
        create_windows=1,
    threads: parameters["threads"]["create_windows"]

    shell:
        " bedtools makewindows -g {input.len} -w {wildcards.window} -s {wildcards.step} > {output.bed} 2>{log.makewin}; "

rule create_bedgraph_track: #
    input:
        track_bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.{track_type}.track.bed",
        windows_bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.win{window}.step{step}.windows.bed",
    output:
        bedgraph=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.{track_type, ^[./]+}.win{window}.step{step}.track.bedgraph"
    log:
        intersect=output_dict["log"]  / "create_bedgraph_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.intersect.log",
        awk=output_dict["log"]  / "create_bedgraph_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.awk.log",
        map=output_dict["log"]  / "create_bedgraph_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.map.log",
        cluster_log=output_dict["cluster_log"] / "create_bedgraph_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_bedgraph_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_bedgraph_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["create_bedgraph_track"],
        time=parameters["time"]["create_bedgraph_track"],
        mem=parameters["memory_mb"]["create_bedgraph_track"],
        create_windows=1,
    threads: parameters["threads"]["create_bedgraph_track"]

    shell:
        " bedtools intersect -wao -a {input.windows_bed}  -b {input.track_bed}  2>{log.intersect} | "
        " awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$NF}}' 2>{log.awk} | "
        " bedtools map -c 4 -o sum -a {input.windows_bed} -b stdin > {output.bedgraph} 2>{log.map} "

rule draw_track: #
    input:
        bedgraph=rules.create_bedgraph_track.output.bedgraph,
        whitelist=rules.select_long_scaffolds.output.whitelist,
        orderlist=rules.select_long_scaffolds.output.orderlist,
        len_file=rules.create_curation_input_links.output.len,
    output:
        png=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.{track_type, ^[./]+}.win{window}.step{step}.png"
    log:
        draw=output_dict["log"]  / "draw_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.draw.log",
        cluster_log=output_dict["cluster_log"] / "draw_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "draw_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "draw_track.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.{track_type}.win{window}.step{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        cpus=parameters["threads"]["draw_track"],
        time=parameters["time"]["draw_track"],
        mem=parameters["memory_mb"]["draw_track"],
    threads: parameters["threads"]["draw_track"]

    shell:
        " PREFIX={output.png}; "
        " PREFIX=${{PREFIX%.png}}; "
        " draw_variant_window_densities.py -i {input.bedgraph} -t bedgraph -o ${{PREFIX}} -l {wildcards.track_type} "
        " -w {wildcards.window} -s {wildcards.step} --density_multiplier 1 "
        " -a {input.whitelist} -n {input.len_file} -z {input.orderlist} --density_thresholds 0.0,0.1,0.4,0.7 "
        " --hide_track_label --rounded --subplots_adjust_left 0.15 > {log.draw} 2>&1; "