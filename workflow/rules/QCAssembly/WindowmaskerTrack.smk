localrules: copy_windowmasker_track

rule windowmasker: #
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
    output:
        counts="{fasta_dir}/repeats/{fasta_prefix, [^/]+}/{fasta_prefix}.windowmasker.counts",
        interval="{fasta_dir}/repeats/{fasta_prefix, [^/]+}/{fasta_prefix}.windowmasker.intervals",
        bed="{fasta_dir}/repeats/{fasta_prefix, [^/]+}/{fasta_prefix}.windowmasker.track.bed",
    log:
        std="{fasta_dir}/windowmasker.{fasta_prefix}.windowmasker.stage1.log",
        cluster_log="{fasta_dir}/windowmasker.{fasta_prefix}.windowmasker.cluster.log",
        cluster_err="{fasta_dir}/windowmasker.{fasta_prefix}.windowmasker.cluster.err"
    benchmark:
        "{fasta_dir}/windowmasker.{fasta_prefix}.windowmasker.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
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
        bed="{fasta_dir}/repeats/{fasta_prefix}/{fasta_prefix}.windowmasker.track.bed",
    output:
        qc_track_bed="{fasta_dir}/assembly_qc/windowmasker/{fasta_prefix, [^/]+}/{fasta_prefix}.windowmasker.track.bed"
    log:
        std="{fasta_dir}/copy_windowmasker_track.{fasta_prefix}.windowmasker.log",
        cluster_log="{fasta_dir}/copy_windowmasker_track.{fasta_prefix}.windowmasker.cluster.log",
        cluster_err="{fasta_dir}/copy_windowmasker_track.{fasta_prefix}.windowmasker.cluster.err"
    benchmark:
        "{fasta_dir}/copy_windowmasker_track.{fasta_prefix}.windowmasker.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("windowmasker"),
        cpus=parameters["threads"]["windowmasker"] ,
        time=parameters["time"]["windowmasker"],
        mem=parameters["memory_mb"]["windowmasker"]
    threads: parameters["threads"]["windowmasker"]
    shell:
        " cp {input.bed} {output.qc_track_bed} > {log.std} 2>&1; "