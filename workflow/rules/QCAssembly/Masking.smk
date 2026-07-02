
rule maskfasta: #
    input:
        fasta="{fasta_dir}/{fasta_prefix}.fasta",
        trf_bed=rules.trf.output.bed if not config["skip_trf"] else [], # trf sometimes hangs on specific genomes
        windowmasker_bed=rules.windowmasker.output.bed
    output:
        masked_fasta="{fasta_dir}/repeats/{fasta_prefix}.softmasked.fasta",
        merged_bed="{fasta_dir}/repeats/{fasta_prefix}.repeats.track.bed",
    log:
        std="{fasta_dir}/maskfasta.{fasta_prefix}.log",
        cluster_log="{fasta_dir}/maskfasta.{fasta_prefix}.cluster.log",
        cluster_err="{fasta_dir}/maskfasta.{fasta_prefix}.cluster.err"
    benchmark:
        "{fasta_dir}/maskfasta.{fasta_prefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"]["name"],
        node_options=parse_node_list("maskfasta"),
        cpus=parameters["threads"]["maskfasta"] ,
        time=parameters["time"]["maskfasta"],
        mem=parameters["memory_mb"]["maskfasta"]
    threads: parameters["threads"]["maskfasta"]
    shell:
        " cat {input.trf_bed} {input.windowmasker_bed} 2>{log.std} | "
        " sort -k1,1V -k2,2n -k3,3n 2>>{log.std} | "
        " bedtools merge -i stdin > {output.merged_bed} 2>>{log.std}; "
        " bedtools maskfasta -soft -fi {input.fasta} -bed {output.merged_bed} -fo {output.masked_fasta} >> {log.std} 2>&1; "
