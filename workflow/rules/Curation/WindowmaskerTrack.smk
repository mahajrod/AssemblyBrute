
rule windowmasker: #
    input:
        fasta=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/{seq_type}/{genome_prefix}.input.{haplotype}.fasta"
    output:
        counts=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/{seq_type}/{genome_prefix}.input.{haplotype}.windowmasker.counts",
        interval=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/{seq_type}/{genome_prefix}.input.{haplotype}.windowmasker.intervals",
        bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/{seq_type}/{genome_prefix}.input.{haplotype}.windowmasker.track.bed",
    log:
        stage1=output_dict["log"]  / "windowmasker.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.stage1.log",
        stage2=output_dict["log"]  / "windowmasker.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.stage2.log",
        conversion=output_dict["log"]  / "windowmasker.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.conversion.log",
        cluster_log=output_dict["cluster_log"] / "windowmasker.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "windowmasker.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "windowmasker.{prev_stage_parameters}..{curation_parameters}.{seq_type}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("windowmasker"),
        cpus=parameters["threads"]["windowmasker"] ,
        time=parameters["time"]["windowmasker"],
        mem=parameters["memory_mb"]["windowmasker"]
    threads: parameters["threads"]["windowmasker"]
    shell:
        " windowmasker -mk_counts -in {input.fasta} -out {output.counts} > {log.stage1} 2>&1; "
        " windowmasker -ustat {output.counts} -in {input.fasta} -out {output.interval} -dust true > {log.stage2} 2>&1;"
        " workflow/scripts/repeats/convert_windowmasker_output_to_bed.py -i {output.interval} -o {output.bed} > {log.conversion} 2>&1;"

"""
rule create_repeat_bins: #
    input:
        fai=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.fasta.fai",
        repeat_bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype}/input/{genome_prefix}.input.{haplotype}.repeat.bed",
    output:
        genome=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.genome",
        #genome_sorted=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.sorted.genome",
        bins=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.bin.bed",
        repeat_bedgraph=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.repeat.bedgraph",
        repeat_binned_bdgraph=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.repeat.binned.bedgraph",
        repeat_binned_no_dot_bedgraph=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.repeat.binned.no_dot.bedgraph"
    params:
        bin_size=lambda wildcards: stage_dict["curation"]["parameters"][wildcards.prev_stage_parameters + ".." + wildcards.curation_parameters]["option_set"]["bin_size"]
        #parameters["tool_options"]["curation"][wildcards.curation_parameters]["bin_size"],
    log:
        cut=output_dict["log"]  / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cut.log",
        sort1=output_dict["log"]  / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.sort1.log",
        make_windows=output_dict["log"]  / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.make_windows.log",
        bedtools=output_dict["log"]  / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.bedtools.log",
        awk=output_dict["log"]  / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.awk.log",
        sed1=output_dict["log"]  / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.sed1.log",
        sort2=output_dict["log"]  / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.sort2.log",
        #cat=output_dict["log"]  / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cat.log",

        #sort3=output_dict["log"]  / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.sort3.log",
        sed2=output_dict["log"]  / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.sed2.log",
        sort4=output_dict["log"]  / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.sort4.log",
        #bedGraphToBigWig=output_dict["log"]  / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.bedGraphToBigWig.log",

        bedtools_map=output_dict["log"]  / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.bedtools_map.log",
        cluster_log=output_dict["cluster_log"] / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_repeat_bins.{prev_stage_parameters}..{curation_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        cpus=parameters["threads"]["create_repeat_bins"] ,
        time=parameters["time"]["create_repeat_bins"],
        mem=parameters["memory_mb"]["create_repeat_bins"]
    threads: parameters["threads"]["create_repeat_bins"]
    shell:
        " cut -f1,2 {input.fai} 2>{log.cut} | sort -k1 -V > {output.genome} 2>{log.sort1}; "
        " bedtools makewindows -g {output.genome} -w {params.bin_size} > {output.bins} 2>{log.make_windows}; "
        " bedtools intersect -a {output.bins} -b {input.repeat_bed} 2>{log.bedtools} | "
        " awk '{{print $0\"\\t\"$3-$2}}' 2>{log.awk} | "
        " sed 's/\./0/g' 2>{log.sed1} | sort -k1 -V >  {output.repeat_bedgraph} 2>{log.sort2}; "
        #" cat {output.bins} 2>{log.cat} | "
        " bedtools map -a {output.bins} -b {output.repeat_bedgraph} -c 4 -o sum -g {output.genome} > {output.repeat_binned_bdgraph} 2>{log.bedtools_map}; "
        #" sort -k2,2 -nr {output.genome} > {output.genome_sorted} 2>{log.sort3}; "
        " sed 's/\./0/g' {output.repeat_binned_bdgraph} 2>{log.sed2} | sort -k1,1 -k2,2n  > {output.repeat_binned_no_dot_bedgraph} 2> {log.sort4}; "
        #" bedGraphToBigWig $$self{outdir}/density_nodot.bed $$self{outdir}/sorted.genome $$self{hgname}_repeat_density.bw"
        #ucsc-bedgraphtobigwig

"""