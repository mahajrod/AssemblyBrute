localrules: tidk_download_db
localrules: filter_tidk_telomere_tracks_for_pretext


rule tidk_download_db:
    input:
        log_dir=out_dir_path / "log/"
    output:
        tidk_db_flag=out_dir_path / "flags/TIDK_DB_FLAG"
    params:
        use_existing_envs=config["use_existing_envs"]
    log:
        std=out_dir_path / "tidk_download_db.log",
        cluster_log=out_dir_path / "log/tidk_download_db.cluster.log",
        cluster_err=out_dir_path / "log/tidk_download_db.cluster.err",
    benchmark:
        out_dir_path / "benchmark/tidk_download_db.benchmark.txt",
    conda:
        config["conda"]["tidk"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["tidk"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("tidk_download_db"),
        cpus=parameters["threads"]["tidk_download_db"] ,
        time=parameters["time"]["tidk_download_db"],
        mem=parameters["memory_mb"]["tidk_download_db"],
    threads: parameters["threads"]["tidk_download_db"]

    shell:
        " if [[ '{params.use_existing_envs}' != 'True' ]]; "
        "   then "
        "   tidk build > {log.std} 2>&1; "
        "   fi; "
        " > {output.tidk_db_flag} >> {log.std} 2>&1; "

rule tidk_search:
    input:
        tidk_db_flag = out_dir_path / "flags/TIDK_DB_FLAG",
        log_dir="{fasta_dir}/log/",
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        canonical_top_kmer="{fasta_dir}/telomere/{fasta_prefix}/{fasta_prefix}.canonical.top.kmer",
        non_canonical_top_kmer="{fasta_dir}/telomere/{fasta_prefix}/{fasta_prefix}.non_canonical.top.kmer"
    output:
        canonical_tidk_bedgraph="{fasta_dir}/telomere_tidk/{fasta_prefix}/{fasta_prefix}.canonical_tidk.bedgraph",
        non_canonical_tidk_bedgraph="{fasta_dir}/telomere_tidk/{fasta_prefix}/{fasta_prefix}.non_canonical_tidk.bedgraph",
    params:
        window_size=parameters["tool_options"]["assembly_qc"]["telomere_tidk_search"]["window_size"],
    log:
        std="{fasta_dir}/log/tidk_search.{fasta_prefix}.log",
        cluster_log="{fasta_dir}/log/tidk_search.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/log/tidk_search.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/log/tidk_search.{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["tidk"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["tidk"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("tidk_search"),
        cpus=parameters["threads"]["tidk_search"] ,
        time=parameters["time"]["tidk_search"],
        mem=parameters["memory_mb"]["tidk_search"],
    threads: parameters["threads"]["tidk_search"]

    shell:
        " OUT_DIR=`dirname {output.canonical_tidk_bedgraph}`; "
        " if [ -s {input.canonical_top_kmer} ]; "
        "   then "
        "   CANNONICAL_OUT_PREFIX={wildcards.fasta_prefix}.canonical; "
        "   CANONICAL_TEL_KMER=`head -n 1 {input.canonical_top_kmer}`; "
        "   echo \"Seeking for cannonical telomere motif ${{CANONICAL_TEL_KMER}}...\" > {log.std} 2>&1; "
        "   tidk search -s ${{CANONICAL_TEL_KMER}} -d ${{OUT_DIR}} -w {params.window_size} "
        "       -o ${{CANNONICAL_OUT_PREFIX}} -e bedgraph {input.fasta} >> {log.std} 2>&1; "
        "   mv ${{OUT_DIR}}/${{CANNONICAL_OUT_PREFIX}}_telomeric_repeat_windows.bedgraph {output.canonical_tidk_bedgraph}; "
        "   else"
        "   echo \"No cannonical telomere motif. Skipping...\" >> {log.std} 2>&1; "
        "   > {output.canonical_tidk_bedgraph}; "
        "   fi;"
        " if [ -s {input.non_canonical_top_kmer} ]; "
        "   then "
        "   NON_CANNONICAL_OUT_PREFIX={wildcards.fasta_prefix}.non_canonical; "
        "   NON_CANONICAL_TEL_KMER=`head -n 1 {input.non_canonical_top_kmer}`; "
        "   echo \"Seeking for non cannonical telomere motif ${{NON_CANONICAL_TEL_KMER}}...\" >> {log.std} 2>&1; "
        "   tidk search -s ${{NON_CANONICAL_TEL_KMER}} -d ${{OUT_DIR}} -w {params.window_size} "
        "               -o ${{NON_CANNONICAL_OUT_PREFIX}} -e bedgraph {input.fasta} >> {log.std} 2>&1; "
        "   mv ${{OUT_DIR}}/${{NON_CANNONICAL_OUT_PREFIX}}_telomeric_repeat_windows.bedgraph {output.non_canonical_tidk_bedgraph}; "
        "   else "
        "   echo \"No non cannonical telomere motif. Skipping...\" >> {log.std} 2>&1; "
        "   > {output.non_canonical_tidk_bedgraph}; "
        "   fi;"

rule filter_tidk_telomere_tracks_for_pretext:
    input:
        canonical_telo_bedgraph="{fasta_dir}/telomere_tidk/{fasta_prefix}/{fasta_prefix}.canonical_tidk.bedgraph",
        non_canonical_telo_bedgraph="{fasta_dir}/telomere_tidk/{fasta_prefix}/{fasta_prefix}.non_canonical_tidk.bedgraph",
        log_dir="{fasta_dir}/log/",
    output:
        canonical_telo_all_bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix, [^/]+}/{fasta_prefix}.canonical_tidk.telomere.all.bedgraph",
        non_canonical_telo_all_bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix, [^/]+}/{fasta_prefix}.non_canonical_tidk.telomere.all.bedgraph",
        canonical_telo_bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix, [^/]+}/{fasta_prefix}.canonical_tidk.telomere.pretext.bedgraph",
        non_canonical_telo_bedgraph="{fasta_dir}/assembly_qc/tracks/{fasta_prefix, [^/]+}/{fasta_prefix}.non_canonical_tidk.telomere.pretext.bedgraph",
    params:
        min_monomer_number=parameters["tool_options"]["assembly_qc"]["telomere_tidk_search"]["min_monomer_number"]
    log:
        canonical="{fasta_dir}/log/copy_tidk_telomere_track_for_pretext.{fasta_prefix}.canonical.log",
        non_canonical="{fasta_dir}/log/copy_tidk_telomere_track_for_pretext.{fasta_prefix}.non_canonical.log",
        cluster_log="{fasta_dir}/log/copy_tidk_telomere_track_for_pretext.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/log/copy_tidk_telomere_track_for_pretext.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/log/copy_telomere_track_for_pretext.{fasta_prefix}.benchmark.txt"
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
