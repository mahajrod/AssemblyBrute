
ruleorder: create_gc_track > create_bedgraph_track
rule create_gc_track: #
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        gc_bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.{track_type, gc}.win{window}.step{step}.track.bedgraph",
    log:
        gc="{fasta_dir}/log/create_gc_track.{fasta_prefix}.{track_type}.win{window}.step{step}.gc.log",
        cluster_log="{fasta_dir}/log/create_gc_track.{fasta_prefix}.{track_type}.win{window}.step{step}.cluster.log",
        cluster_err="{fasta_dir}/log/create_gc_track.{fasta_prefix}.{track_type}.win{window}.step{step}.cluster.err"
    benchmark:
        "{fasta_dir}/log/create_gc_track.{fasta_prefix}.{track_type}.win{window}.step{step}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("create_gc_track"),
        cpus=parameters["threads"]["create_gc_track"],
        time=parameters["time"]["create_gc_track"],
        mem=parameters["memory_mb"]["create_gc_track"]
    threads: parameters["threads"]["create_gc_track"]

    shell:
        " workflow/scripts/curation/count_gc_in_windows.py -i {input.fasta} "
        "     -w {wildcards.window} -s {wildcards.step} -o {output.gc_bedgraph} 2>{log.gc}; "
