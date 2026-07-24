localrules: copy_windowmasker_track

rule copy_windowmasker_track: #
    input:
        bed="{fasta_dir}/{fasta_prefix}/repeats/{fasta_prefix}.windowmasker.track.bed",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        qc_track_bed="{fasta_dir}/assembly_qc/windowmasker/{fasta_prefix}/{fasta_prefix}.windowmasker.track.bed"
    log:
        std="{fasta_dir}/log/copy_windowmasker_track.{fasta_prefix}.windowmasker.log",
        cluster_log="{fasta_dir}/log/copy_windowmasker_track.{fasta_prefix}.windowmasker.cluster.log",
        cluster_err="{fasta_dir}/log/copy_windowmasker_track.{fasta_prefix}.windowmasker.cluster.err"
    benchmark:
        "{fasta_dir}/log/copy_windowmasker_track.{fasta_prefix}.windowmasker.benchmark.txt"
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