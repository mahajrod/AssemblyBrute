rule blast2_2_27:
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta"
    output:
        nhr="{fasta_dir}/blast2.2.27/{fasta_prefix}.fasta.nhr",
        nin="{fasta_dir}/blast2.2.27/{fasta_prefix}.fasta.nin",
        nsq="{fasta_dir}/blast2.2.27/{fasta_prefix}.fasta.nsq",
    log:
        std="{fasta_dir}/blast2.2.27.{fasta_prefix}.log",
        cluster_log="{fasta_dir}/blast2.2.27.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/blast2.2.27.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/blast2.2.27.{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["simulate_pcr"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["simulate_pcr"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("blast"),
        cpus=parameters["threads"]["blast"] ,
        time=parameters["time"]["blast"],
        mem=parameters["memory_mb"]["blast"]
    threads: parameters["threads"]["blast"]

    shell:
        " makeblastdb -dbtype nucl -in {input.fasta} "
        "             -out {wildcards.fasta_dir}/blast2.2.27/{wildcards.fasta_prefix}.fasta > {log.std} 2>&1;"