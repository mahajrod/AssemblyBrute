localrules: tidk_download_db


rule tidk_download_db:
    input:
        log_dir=ancient(config["out_dir"] / "log/")
    output:
        tidk_db_flag=config["out_dir"] / "flags/TIDK_DB_FLAG"
    params:
        use_existing_envs=config["use_existing_envs"]
    log:
        std=config["out_dir"] / "tidk_download_db.log",
        cluster_log=config["out_dir"] / "log/tidk_download_db.cluster.log",
        cluster_err=config["out_dir"] / "log/tidk_download_db.cluster.err",
    benchmark:
        config["out_dir"] / "benchmark/tidk_download_db.benchmark.txt",
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
        " then "
        "     tidk build > {log.std} 2>&1; "
        " fi; "
        " > {output.tidk_db_flag} >> {log.std} 2>&1; "

rule tidk_search:
    input:
        log_dir=ancient("{fasta_dir}/log/"),
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        canonical_top_kmer="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.canonical.top.kmer",
        non_canonical_top_kmer="{fasta_dir}/{fasta_prefix}/telomere/{fasta_prefix}.non_canonical.top.kmer"
    output:
        canonical_tidk_bedgraph="{fasta_dir}/{fasta_prefix}/telomere_tidk/{fasta_prefix}.canonical_tidk.win{window}.step{window}.bedgraph",
        non_canonical_tidk_bedgraph="{fasta_dir}/{fasta_prefix}/telomere_tidk/{fasta_prefix}.non_canonical_tidk.win{window}.step{window}.bedgraph",
    log:
        std="{fasta_dir}/log/tidk_search.{fasta_prefix}.win{window}.step{window}.log",
        cluster_log="{fasta_dir}/log/tidk_search.{fasta_prefix}.win{window}.step{window}.cluster.log",
        cluster_err="{fasta_dir}/log/tidk_search.{fasta_prefix}.win{window}.step{window}.cluster.err"
    benchmark:
        "{fasta_dir}/log/tidk_search.{fasta_prefix}.win{window}.step{window}.benchmark.txt"
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
        " then "
        "     CANNONICAL_OUT_PREFIX={wildcards.fasta_prefix}.canonical; "
        "     CANONICAL_TEL_KMER=`head -n 1 {input.canonical_top_kmer}`; "
        "     echo \"Seeking for cannonical telomere motif ${{CANONICAL_TEL_KMER}}...\" > {log.std} 2>&1; "
        "     tidk search -s ${{CANONICAL_TEL_KMER}} -d ${{OUT_DIR}} -w {wildcards.window} "
        "         -o ${{CANNONICAL_OUT_PREFIX}} -e bedgraph {input.fasta} >> {log.std} 2>&1; "
        "     mv ${{OUT_DIR}}/${{CANNONICAL_OUT_PREFIX}}_telomeric_repeat_windows.bedgraph {output.canonical_tidk_bedgraph}; "
        " else"
        "     echo \"No cannonical telomere motif. Skipping...\" >> {log.std} 2>&1; "
        "     > {output.canonical_tidk_bedgraph}; "
        " fi;"
        " if [ -s {input.non_canonical_top_kmer} ]; "
        " then "
        "     NON_CANNONICAL_OUT_PREFIX={wildcards.fasta_prefix}.non_canonical; "
        "     NON_CANONICAL_TEL_KMER=`head -n 1 {input.non_canonical_top_kmer}`; "
        "     echo \"Seeking for non cannonical telomere motif ${{NON_CANONICAL_TEL_KMER}}...\" >> {log.std} 2>&1; "
        "     tidk search -s ${{NON_CANONICAL_TEL_KMER}} -d ${{OUT_DIR}} -w {wildcards.window} "
        "         -o ${{NON_CANNONICAL_OUT_PREFIX}} -e bedgraph {input.fasta} >> {log.std} 2>&1; "
        "     mv ${{OUT_DIR}}/${{NON_CANNONICAL_OUT_PREFIX}}_telomeric_repeat_windows.bedgraph {output.non_canonical_tidk_bedgraph}; "
        " else "
        "     echo \"No non cannonical telomere motif. Skipping...\" >> {log.std} 2>&1; "
        "     > {output.non_canonical_tidk_bedgraph}; "
        " fi;"
