rule quast:
    input:
        assembly="{fasta_dir}/{fasta_prefix}.fasta",
        log_dir=ancient("{fasta_dir}/log/"),
    output:
        dir=directory("{fasta_dir}/assembly_qc/quast/{fasta_prefix}"),
    params:
        large_genome_flag="--large" if parameters["tool_options"]["quast"]["large_genome"] else "",
    log:
        std="{fasta_dir}/log/quast.{fasta_prefix}.log",
        cluster_log="{fasta_dir}/log/quast.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/log/quast.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/log/quast.{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("quast"),
        cpus=parameters["threads"]["quast"],
        time=parameters["time"]["quast"],
        mem=parameters["memory_mb"]["quast"],
    threads:
        parameters["threads"]["quast"]
    shell:
         " quast -o {output.dir} -t {threads} -l {wildcards.fasta_prefix} {params.large_genome_flag} "
         "     {input.assembly} > {log.std} 2>&1;"
