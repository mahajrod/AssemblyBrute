
localrules: filter_tidk_telomere_tracks_for_pretext

rule filter_tidk_telomere_tracks_for_pretext:
    input:
        canonical_telo_bedgraph="{fasta_dir}/{fasta_prefix}/telomere_tidk/{fasta_prefix}.canonical_tidk.win{window}.step{window}.bedgraph",
        non_canonical_telo_bedgraph="{fasta_dir}/{fasta_prefix}/telomere_tidk/{fasta_prefix}.non_canonical_tidk.win{window}.step{window}.bedgraph",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        canonical_telo_all_bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.canonical_tidk_telomere_all.win{window}.step{window}.track.bedgraph",
        non_canonical_telo_all_bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.non_canonical_tidk_telomere_all.win{window}.step{window}.track.bedgraph",
        canonical_telo_bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.canonical_tidk_telomere_filtered.win{window}.step{window}.track.bedgraph",
        non_canonical_telo_bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix}/{fasta_prefix}.non_canonical_tidk_telomere_filtered.win{window}.step{window}.track.bedgraph",
    params:
        min_monomer_number=parameters["tool_options"]["assembly_qc"]["telomere_tidk_search"]["min_monomer_number"]
    log:
        canonical="{fasta_dir}/log/copy_tidk_telomere_track_for_pretext.{fasta_prefix}.win{window}.step{window}.canonical.log",
        non_canonical="{fasta_dir}/log/copy_tidk_telomere_track_for_pretext.{fasta_prefix}.win{window}.step{window}.non_canonical.log",
        cluster_log="{fasta_dir}/log/copy_tidk_telomere_track_for_pretext.{fasta_prefix}.win{window}.step{window}.cluster.log",
        cluster_err="{fasta_dir}/log/copy_tidk_telomere_track_for_pretext.{fasta_prefix}.win{window}.step{window}.cluster.err"
    benchmark:
        "{fasta_dir}/log/copy_telomere_track_for_pretext.{fasta_prefix}.win{window}.step{window}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("telo_container"),
        cpus=parameters["threads"]["copy_telomere_track_for_pretext"] ,
        time=parameters["time"]["copy_telomere_track_for_pretext"],
        mem=parameters["memory_mb"]["copy_telomere_track_for_pretext"],
    threads: parameters["threads"]["copy_telomere_track_for_pretext"]

    shell:
        " cp -f {input.canonical_telo_bedgraph} {output.canonical_telo_all_bedgraph} > {log.canonical} 2>&1; "
        " cp -f {input.non_canonical_telo_bedgraph} {output.non_canonical_telo_all_bedgraph} > {log.non_canonical} 2>&1; "
        " awk '{{ if ($4 >= {params.min_monomer_number}) print $0}}' {input.canonical_telo_bedgraph} > {output.canonical_telo_bedgraph} 2>> {log.canonical}; "
        " awk '{{ if ($4 >= {params.min_monomer_number}) print $0}}' {input.non_canonical_telo_bedgraph} > {output.non_canonical_telo_bedgraph} 2>> {log.non_canonical}; "
