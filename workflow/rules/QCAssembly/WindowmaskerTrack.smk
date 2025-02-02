localrules: copy_windowmasker_track

rule windowmasker: #
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
    output:
        counts="{fasta_dir}/repeats/{fasta_prefix, [^/]+}.{track_type, windowmasker}.counts",
        interval="{fasta_dir}/repeats/{fasta_prefix, [^/]+}.{track_type, windowmasker}.intervals",
        bed="{fasta_dir}/repeats/{fasta_prefix, [^/]+}.{track_type, windowmasker}.track.bed",
        #qc_track_bed="{fasta_dir}/assembly_qc/{track_type, windowmasker}/{fasta_prefix, [^/]+}/{fasta_prefix}.{track_type}.track.bed"
    log:
        std="{fasta_dir}/windowmasker.{fasta_prefix}.{track_type}.stage1.log",
        cluster_log="{fasta_dir}/windowmasker.{fasta_prefix}.{track_type}.cluster.log",
        cluster_err="{fasta_dir}/windowmasker.{fasta_prefix}.{track_type}.cluster.err"
    benchmark:
        "{fasta_dir}/windowmasker.{fasta_prefix}.{track_type}.benchmark.txt"
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
        " windowmasker -mk_counts -in {input.fasta} -out {output.counts} > {log.std} 2>&1; "
        " windowmasker -ustat {output.counts} -in {input.fasta} -out {output.interval} -dust true >> {log.std} 2>&1; "
        " workflow/scripts/repeats/convert_windowmasker_output_to_bed.py -i {output.interval} -o {output.bed} >> {log.std} 2>&1; "


rule copy_windowmasker_track: #
    input:
        bed="{fasta_dir}/repeats/{fasta_prefix}.{track_type}.track.bed",
    output:
        qc_track_bed="{fasta_dir}/assembly_qc/{track_type, windowmasker}/{fasta_prefix, [^/]+}/{fasta_prefix}.{track_type}.track.bed"
    log:
        std="{fasta_dir}/copy_windowmasker_track.{fasta_prefix}.{track_type}log",
        cluster_log="{fasta_dir}/copy_windowmasker_track.{fasta_prefix}.{track_type}.cluster.log",
        cluster_err="{fasta_dir}/copy_windowmasker_track.{fasta_prefix}.{track_type}.cluster.err"
    benchmark:
        "{fasta_dir}/copy_windowmasker_track.{fasta_prefix}.{track_type}.benchmark.txt"
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