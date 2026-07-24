
rule windowmasker: #
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        counts="{fasta_dir}/{fasta_prefix}/repeats/{fasta_prefix}.windowmasker.counts",
        interval="{fasta_dir}/{fasta_prefix}/repeats/{fasta_prefix}.windowmasker.intervals",
        bed="{fasta_dir}/{fasta_prefix}/repeats/{fasta_prefix}.windowmasker.track.bed",
    log:
        std="{fasta_dir}/log/windowmasker.{fasta_prefix}.windowmasker.stage1.log",
        cluster_log="{fasta_dir}/log/windowmasker.{fasta_prefix}.windowmasker.cluster.log",
        cluster_err="{fasta_dir}/log/windowmasker.{fasta_prefix}.windowmasker.cluster.err"
    benchmark:
        "{fasta_dir}/log/windowmasker.{fasta_prefix}.windowmasker.benchmark.txt"
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
