
rule windowmasker: #
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
    output:
        counts="{fasta_dir}/repeats/{fasta_prefix}.{track_type, windowmasker}.counts",
        interval="{fasta_dir}/repeats/{fasta_prefix}.{track_type, windowmasker}.intervals",
        bed="{fasta_dir}/repeats/{fasta_prefix}.{track_type, windowmasker}.track.bed",
        qc_track_bed="{fasta_dir}/assembly_qc/{track_type, trf}/{fasta_prefix, [^/]+}/{fasta_prefix}.{track_type, windowmasker}.track.bed"
    log:
        std="{fasta_dir}/{fasta_prefix}.{track_type}.stage1.log",
        cluster_log="{fasta_dir}/{fasta_prefix}.{track_type}.cluster.log",
        cluster_err="{fasta_dir}/{fasta_prefix}.{track_type}.cluster.err"
    benchmark:
        "{fasta_dir}/{fasta_prefix}.{track_type}.benchmark.txt"
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
        " LOG=`realpath {log.std}`; "
        " windowmasker -mk_counts -in {input.fasta} -out {output.counts} > ${{LOG}} 2>&1; "
        " windowmasker -ustat {output.counts} -in {input.fasta} -out {output.interval} -dust true > ${{LOG}} 2>&1; "
        " workflow/scripts/repeats/convert_windowmasker_output_to_bed.py -i {output.interval} -o {output.bed} > ${{LOG}} 2>&1; "
        " cp {output.bed} {output.qc_track_bed} > ${{LOG}} 2>&1; "
