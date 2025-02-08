
rule trf: #
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
    output:
        simple_bed="{fasta_dir}/repeats/{fasta_prefix}/{fasta_prefix, [^/]+}.trf.simple.bed",
        bed="{fasta_dir}/repeats/{fasta_prefix}/{fasta_prefix, [^/]+}.trf.track.bed",
        #qc_track_bed="{fasta_dir}/assembly_qc/{track_type, trf}/{fasta_prefix, [^/]+}/{fasta_prefix}.{track_type, trf}.track.bed"
    params:
        matching_weight=parse_option("matching_weight", parameters["tool_options"]["trf"], " -m "),
        mismatching_penalty=parse_option("mismatching_penalty", parameters["tool_options"]["trf"], " -s "),
        indel_penalty=parse_option("indel_penalty", parameters["tool_options"]["trf"], " -l "),
        match_probability=parse_option("match_probability", parameters["tool_options"]["trf"], " -a "),
        indel_probability=parse_option("indel_probability", parameters["tool_options"]["trf"], " -d "),
        min_alignment_score=parse_option("min_alignment_score", parameters["tool_options"]["trf"], " -c "),
        max_period=parse_option("max_period", parameters["tool_options"]["trf"], " -e "),
        max_repeat_length=parse_option("max_repeat_length", parameters["tool_options"]["trf"], " -g "),
    log:
        std="{fasta_dir}/trf.{fasta_prefix}.trf.log",
        cluster_log="{fasta_dir}/trf.{fasta_prefix}.trf.cluster.log",
        cluster_err="{fasta_dir}/trf.{fasta_prefix}.trf.cluster.err"
    benchmark:
        "{fasta_dir}/trf.{fasta_prefix}.trf.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("trf"),
        cpus=parameters["threads"]["trf"] ,
        time=parameters["time"]["trf"],
        mem=parameters["memory_mb"]["trf"]
    threads: parameters["threads"]["trf"]
    shell:
        " LOG=`realpath {log.std}`; "
        " WORK_DIR=`dirname {output.bed}`; "
        " INPUT_FASTA=`basename  {input.fasta}`; "
        " OUTPUT_PREFIX=`basename {output.simple_bed}`; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.simple.bed}};"
        " cd ${{WORK_DIR}}; "
        " tandem_repeat_masking.py --sleep 120 -t {threads} {params.matching_weight} {params.mismatching_penalty} "
        " {params.indel_penalty} {params.match_probability} {params.indel_probability} {params.min_alignment_score} "
        " {params.max_period} {params.max_repeat_length} -i ../${{INPUT_FASTA}} -o ${{OUTPUT_PREFIX}} > ${{LOG}} 2>&1; "
        " cut -f1-3 `basename {output.simple_bed}` 2>>${{LOG}} | grep -vP '^#' 2>>${{LOG}} | sort -k1,1V -k2,2n -k3,3n  2>>${{LOG}} | "
        " bedtools merge -i stdin > `basename {output.bed}` 2>>${{LOG}}; "


rule copy_trf_track: #
    input:
        bed="{fasta_dir}/repeats/{fasta_prefix}/{fasta_prefix}.trf.track.bed",
    output:
        qc_track_bed="{fasta_dir}/{fasta_prefix}/assembly_qc/trf/{fasta_prefix, [^/]+}/{fasta_prefix}.trf.track.bed"
    log:
        std="{fasta_dir}/copy_trf_track.{fasta_prefix}.trf.log",
        cluster_log="{fasta_dir}/copy_trf_track.{fasta_prefix}.trf.cluster.log",
        cluster_err="{fasta_dir}/copy_trf_track.{fasta_prefix}.trf.cluster.err"
    benchmark:
        "{fasta_dir}/copy_trf_track.{fasta_prefix}.trf.benchmark.txt"
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
        " cp {input.bed} {output.qc_track_bed} > {log.std} 2>&1; "
