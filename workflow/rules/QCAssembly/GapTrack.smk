

rule create_gap_track: #
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        gap_bed="{fasta_dir}/assembly_qc/{track_type, gap}/{fasta_prefix}/{fasta_prefix}.{track_type}.track.bed",
        gap_bedgraph="{fasta_dir}/assembly_qc/{track_type, gap}/{fasta_prefix}/{fasta_prefix}.{track_type}.track.bedgraph",
        gap_bedgraph_alias="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.{track_type, gap}.track.bedgraph",
    log:
        seqtk="{fasta_dir}/log/create_gap_track.{fasta_prefix}.{track_type}.seqtk.log",
        awk="{fasta_dir}/log/create_gap_track.{fasta_prefix}.{track_type}.awk.log",
        cp="{fasta_dir}/log/create_gap_track.{fasta_prefix}.{track_type}.cp.log",
        cluster_log="{fasta_dir}/log/create_gap_track.{fasta_prefix}.{track_type}.cluster.log",
        cluster_err="{fasta_dir}/log/create_gap_track.{fasta_prefix}.{track_type}.cluster.err"
    benchmark:
        "{fasta_dir}/log/create_gap_track.{fasta_prefix}.{track_type}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_gap_track"),
        cpus=parameters["threads"]["create_gap_track"],
        time=parameters["time"]["create_gap_track"],
        mem=parameters["memory_mb"]["create_gap_track"]
    threads: parameters["threads"]["create_gap_track"]

    shell:
        " seqtk cutN -n 1 -g  {input.fasta} > {output.gap_bed} 2>{log.seqtk}; "
        " awk '{{printf \"%s\\t%i\\n\",$0,1}}' {output.gap_bed} > {output.gap_bedgraph} 2>{log.awk}; "
        " cp {output.gap_bedgraph} {output.gap_bedgraph_alias} > {log.cp} 2>&1; "
