

rule ragtag: #
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        reference_fasta=config["out_dir"] / "data/reference/{reference}/{reference}/masking/{reference}.softmasked.fasta",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        ragtag_fasta="{fasta_dir}/{fasta_prefix}/ragtag/{reference}/{fasta_prefix}.to.{reference}.fasta",
        alias_ragtag_fasta="{fasta_dir}/{fasta_prefix}/ragtag/{reference}/ragtag.scaffold.fasta",
        ragtag_agp="{fasta_dir}/{fasta_prefix}/ragtag/{reference}/{fasta_prefix}.to.{reference}.agp",
        alias_ragtag_agp="{fasta_dir}/{fasta_prefix}/ragtag/{reference}/ragtag.scaffold.agp",
        ragtag_stats="{fasta_dir}/{fasta_prefix}/ragtag/{reference}/{fasta_prefix}.to.{reference}.stats",
        alias_stats="{fasta_dir}/{fasta_prefix}/ragtag/{reference}/ragtag.scaffold.stats",
    log:
        ragtag="{fasta_dir}/log/ragtag.{fasta_prefix}.{reference}.ragtag.log",
        ln="{fasta_dir}/log/ragtag.{fasta_prefix}.{reference}.ln.log",
        cluster_log="{fasta_dir}/log/ragtag.{fasta_prefix}.{reference}.cluster.log",
        cluster_err="{fasta_dir}/log/ragtag.{fasta_prefix}.{reference}.cluster.err"
    benchmark:
        "{fasta_dir}/log/ragtag.{fasta_prefix}.{reference}.benchmark.txt"
    conda:
        config["conda"]["ragtag"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["ragtag"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("ragtag"),
        cpus=parameters["threads"]["ragtag"],
        time=parameters["time"]["ragtag"],
        mem=parameters["memory_mb"]["ragtag"],
        ragtag=1
    threads: parameters["threads"]["ragtag"]

    shell:
        " RAGTAG_DIR=`dirname {output.ragtag_fasta}`; "
        " ragtag.py scaffold -t {threads} -o ${{RAGTAG_DIR}}  -w {input.reference_fasta} {input.fasta} > {log.ragtag} 2>&1; "
        " ln -sf ragtag.scaffold.fasta {output.ragtag_fasta} > {log.ln} 2>&1; "
        " ln -sf ragtag.scaffold.agp {output.ragtag_agp} > {log.ln} 2>&1; "
        " ln -sf ragtag.scaffold.stats {output.ragtag_stats} > {log.ln} 2>&1; "